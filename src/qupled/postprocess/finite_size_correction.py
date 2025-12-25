from __future__ import annotations
import math

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from qupled.database.database_handler import DataBaseHandler
from qupled.database.finite_size_correction_tables import FiniteSizeCorrectionTables
from qupled.database.scheme_tables import BaseTableKeys, TableKeys, RunStatus
from qupled.schemes import hf
from qupled.util.dimension import Dimension
from qupled.util import serialize


class FiniteSizeCorrectionError(Exception):
    """Custom exception for finite size correction operations."""

    pass


FLOAT_TOLERANCE = 1e-12


class Run:
    def __init__(self, rs, id):
        self.rs: float = rs
        self.id: int = id


class FiniteSizeCorrection:
    """
    Post-processing FSC solver. It computes and writes results
    into the *scheme* run at (rs=target_rs, theta, dimension).
    """

    def __init__(self):
        self.db_handler = DataBaseHandler()
        self.results: Result = Result()
        self.solver: hf.Solver | None = None
        self.inputs: Input | None = None
        self.runs: dict[int, Run] = {}
        self.target_run_id: int | None = None
        self.status: RunStatus | None = None

    @property
    def run_id(self):
        """
        Property that retrieves the run ID from the database handler.

        Returns:
            str: The run ID associated with the current database handler.
        """
        return self._db_tables.run_id if self._db_tables is not None else None

    def compute(self, solver: hf.Solver, inputs: Input):
        """
        Compute the finite size correction and store the results in the database.

        Args:
            inputs: Input parameters.
        """
        self.solver = solver
        self.inputs = inputs
        self._add_run_to_database()
        self._compute_correction()
        self._save()

    @property
    def _db_tables(self) -> FiniteSizeCorrectionTables | None:
        """
        Retrieves the finite size correction tables associated with the current database handler.
        Returns:
            FiniteSizeCorrectionTables: The finite size correction tables from the database handler.
        """
        return self.db_handler.fsc_tables if self.db_handler is not None else None

    def _compute_correction(self):
        """
        Computes the finite size correction by solving the scheme for all coupling values,
        calculating the correction, and updating the results and status.
        """
        try:
            self._solve_scheme_for_all_coupling()
            correction = Correction(self.db_handler)
            correction.compute(self.runs, self.inputs)
            self.results = Result(
                uint=correction.internal_energy,
                fxc=correction.free_energy,
            )
            self.status = RunStatus.SUCCESS
        except FiniteSizeCorrectionError as e:
            print(f"Error during finite size correction computation: {e}")
            self.status = RunStatus.FAILED

    def _add_run_to_database(self):
        """
        Adds the current run to the database using the input parameters.
        """
        self._db_tables.insert_run(self.inputs)

    def _save(self):
        """
        Saves the current state and results to the database. Updates the run status,
        scheme run ID, and inserts the computed results.
        """
        self._db_tables.update_run_status(self.status)
        self._db_tables.update_scheme_run_id(self.target_run_id)
        self._db_tables.insert_results(self.results.__dict__)

    def _solve_scheme_for_all_coupling(self):
        """
        Solves the scheme for all coupling values in the grid. Updates the runs dictionary
        and sets the target run ID.
        """
        scheme_input = self.inputs.scheme
        rs_grid = self._build_rs_grid()
        self._search_results_in_the_database(rs_grid)
        filtered_rs_grid = [
            rs for rs in rs_grid if self._key_from_float(rs) not in self.runs
        ]
        target_coupling = scheme_input.coupling
        for rs in filtered_rs_grid:
            self._solve_scheme_for_one_coupling(rs)
        scheme_input.coupling = target_coupling
        target_key = self._key_from_float(target_coupling)
        self.target_run_id = self.runs[target_key].id
        print(f"--- Runs completed")

    def _build_rs_grid(self) -> np.ndarray:
        """
        Builds the grid of coupling parameter (rs) values based on the input resolution.

        Returns:
            np.ndarray: Array of coupling parameter values.
        """
        drs = self.inputs.drs
        if drs <= 0:
            raise FiniteSizeCorrectionError(
                f"Invalid coupling resolution, it must be at greater than 0"
            )
        target = self.inputs.scheme.coupling
        grid = np.arange(0.0, target + drs, drs)
        if len(grid) == 0 or not math.isclose(
            grid[-1], target, rel_tol=0, abs_tol=FLOAT_TOLERANCE
        ):
            grid = np.append(grid, target)
        if not math.isclose(grid[0], 0.0, rel_tol=0, abs_tol=FLOAT_TOLERANCE):
            grid = np.insert(grid, 0, 0.0)
        return grid.astype(float)

    def _solve_scheme_for_one_coupling(self, rs: float) -> None:
        """
        Solves the scheme for a single coupling value.

        Args:
            rs (float): Coupling parameter value.
        """
        print(f"--- Computing run for rs={rs}")
        self.inputs.scheme.coupling = rs
        self.solver.compute(self.inputs.scheme)
        run_status = self.solver.get_solver_status()
        if run_status != RunStatus.SUCCESS:
            raise FiniteSizeCorrectionError(
                f"Run for rs={rs} failed with status {run_status}"
            )
        key = self._key_from_float(rs)
        self.runs[key] = Run(rs, self.solver.run_id)

    def _key_from_float(self, x):
        """
        Converts a float value to an integer key based on the input resolution.

        Args:
            x (float): Float value to convert.

        Returns:
            int: Integer key.
        """
        precision = int(np.abs(np.log10(self.inputs.drs)))
        return int(round(x * 10**precision))

    def _search_results_in_the_database(self, rs_grid):
        """
        Searches for existing results in the database for the given coupling grid.

        Args:
            rs_grid (list[float]): List of coupling parameter values.
        """
        inputs_scheme = self.inputs.scheme
        scheme_tables = self.db_handler.scheme_tables
        runs = scheme_tables.inspect_runs()
        filtered_runs = [
            run
            for run in runs
            if run[TableKeys.THEORY.value] == inputs_scheme.theory
            and run[TableKeys.DEGENERACY.value] == inputs_scheme.degeneracy
        ]
        for rs in rs_grid:
            # It's better not to use floats as keys in a dictionary
            key = self._key_from_float(rs)
            run = next(
                (
                    run
                    for run in filtered_runs
                    if self._key_from_float(run[TableKeys.COUPLING.value]) == key
                ),
                None,
            )
            if run is None:
                continue
            run_id = run[BaseTableKeys.PRIMARY_KEY.value]
            run_inputs = scheme_tables.get_inputs(run_id)
            if (
                run_inputs["cutoff"] == inputs_scheme.cutoff
                and run_inputs["matsubara"] == inputs_scheme.matsubara
                and run_inputs["resolution"] == inputs_scheme.resolution
            ):
                self.runs[key] = Run(rs, run_id)


@serialize.serializable_dataclass
class Input:
    """
    Class used to store the inputs for the :obj:`qupled.finite_size_correction.FiniteSizeCorrection` class.
    """

    number_of_particles: int
    """Number of particles."""
    drs: float = 0.1
    """Coupling parameter resolution. Default = ``0.1``"""
    scheme: hf.Input = None
    """Scheme inputs. Default = ``None``"""


@serialize.serializable_dataclass
class Result:
    uint: float | None = None
    """Internal energy"""
    fxc: float | None = None
    """Free energy"""


class Correction:
    def __init__(self, db_handler: DataBaseHandler):
        """
        Initialize a Correction instance.

        Args:
            db_handler (DataBaseHandler): Database handler for accessing data.
        """
        self.continuous: dict[int, list[float]] = {}
        self.db_handler = db_handler
        self.free_energy: float | None = None
        self.inputs: Input | None = None
        self.internal_energy: float | None = None
        self.rs_grid: list[float] = []
        self.runs: dict[int, Run] = {}
        self.ssfi: dict[int, interp1d] = {}

    def compute(self, runs: dict[int, Run], inputs: Input):
        """
        Computes the finite size correction using the provided runs and inputs.

        Args:
            runs (dict[int, Run]): Dictionary of runs with their coupling values.
            inputs (Input): Input parameters for the correction computation.
        """
        self.runs = runs
        self.inputs = inputs
        self._build_rs_grid()
        self._build_interpolators()
        self._compute_continous()
        self._compute_discrete()

    def _build_rs_grid(self):
        """
        Builds the grid of coupling parameter (rs) values from the runs.
        """
        if not self.rs_grid:
            for _, run in self.runs.items():
                self.rs_grid.append(run.rs)
                self.rs_grid.sort()

    def _build_interpolators(self):
        """
        Builds interpolators for the structure factor data from the runs.
        """
        if not self.ssfi:
            for key, run in self.runs.items():
                wvg, ssf = self._read_wvg_ssf(run.id)
                self.ssfi[key] = interp1d(
                    wvg, ssf, kind="cubic", fill_value="extrapolate"
                )

    def _read_wvg_ssf(self, run_id: int) -> tuple[list[float], list[float]]:
        """
        Reads the wavevector (wvg) and structure factor (ssf) data for a given run ID.

        Args:
            run_id (int): Run ID to fetch data for.

        Returns:
            tuple[list[float], list[float]]: Wavevector and structure factor data.
        """
        data = self.db_handler.scheme_tables.get_run(run_id)
        res = data.get("results", {})
        wvg = res.get("wvg", None)
        ssf = res.get("ssf", None)
        if wvg is None or ssf is None:
            raise FiniteSizeCorrectionError(f"Malformed results for Run {run_id}.")
        return wvg, ssf

    def _is_2D(self):
        """
        Checks if the current scheme is 2D.

        Returns:
            bool: True if the scheme is 2D, False otherwise.
        """
        return self.inputs.scheme.dimension == Dimension._2D

    def _compute_continous(self):
        """
        Computes the continuous part of the finite size correction.
        """
        inputs_scheme = self.inputs.scheme
        if not self.continuous:
            cutoff = inputs_scheme.cutoff
            for key, _ in self.runs.items():
                ssfi = self.ssfi[key]
                self.continuous[key] = (
                    self._compute_continuous_2D(ssfi, cutoff)
                    if self._is_2D()
                    else self._compute_continuous_3D(ssfi, cutoff)
                )

    def _compute_continuous_3D(self, ssfi: interp1d, cutoff: float) -> float:
        """
        Computes the continuous part of the correction for a 3D system.

        Args:
            ssfi (interp1d): Interpolator for the structure factor.
            cutoff (float): Cutoff value for the integral.

        Returns:
            float: Computed continuous correction value.
        """
        _lambda = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
        val = quad(lambda q: float(ssfi(q) - 1.0), 0.0, cutoff, limit=200)[0]
        return val / (np.pi * _lambda)

    def _compute_continuous_2D(self, ssfi: interp1d, cutoff: float) -> float:
        """
        Computes the continuous part of the correction for a 2D system.

        Args:
            ssfi (interp1d): Interpolator for the structure factor.
            cutoff (float): Cutoff value for the integral.

        Returns:
            float: Computed continuous correction value.
        """
        _lambda = 1.0 / np.sqrt(2.0)
        val = quad(lambda q: float(ssfi(q) - 1.0), 0.0, cutoff, limit=200)[0]
        return val / (2.0 * _lambda)

    def _compute_discrete(self):
        """
        Computes the discrete part of the finite size correction.
        """
        rs_uint_vals = []
        rs_uint_target: float | None = None
        inputs_scheme = self.inputs.scheme
        number_of_particles = self.inputs.number_of_particles
        cutoff = self.inputs.scheme.cutoff
        for key, run in self.runs.items():
            ssfi = self.ssfi[key]
            discrete = (
                self._compute_discrete_2D(ssfi, cutoff, number_of_particles)
                if self._is_2D()
                else self._compute_discrete_3D(ssfi, cutoff, number_of_particles)
            )
            madelung = self._compute_madelung(number_of_particles)
            rs_uint = self.continuous[key] - discrete - 0.5 * madelung
            rs_uint_vals.append(rs_uint)
            if math.isclose(
                run.rs,
                inputs_scheme.coupling,
                rel_tol=0,
                abs_tol=FLOAT_TOLERANCE,
            ):
                rs_uint_target = rs_uint
        if rs_uint_target is None:
            raise FiniteSizeCorrectionError(
                "Target coupling parameter not found in coupling grid"
            )
        self._compute_output(rs_uint_vals, rs_uint_target)

    def _compute_output(self, rs_uint, rs_uint_target):
        """
        Computes the final output values for internal and free energy corrections.

        Args:
            rs_uint (list[float]): List of internal energy values for each coupling.
            rs_uint_target (float): Internal energy value at the target coupling.
        """
        interp_kind = "cubic" if len(self.runs) >= 4 else "linear"
        coupling = self.inputs.scheme.coupling
        # uint FSC at target
        self.internal_energy = rs_uint_target / coupling
        # fxc FSC: integrate rs*uint(rs) from 0 to target_rs, then divide by rs^2
        rs_grid_np = np.asarray(self.rs_grid, dtype=float)
        rs_uint_np = np.asarray(rs_uint, dtype=float)
        interp = interp1d(
            rs_grid_np, rs_uint_np, kind=interp_kind, fill_value="extrapolate"
        )
        integral_val = quad(lambda r: float(interp(r)), 0.0, coupling, limit=200)[0]
        self.free_energy = integral_val / (coupling**2.0)

    def _compute_discrete_2D(
        ssfi: interp1d, qcap: float, number_of_particles: int
    ) -> float:
        """
        Computes the discrete part of the correction for a 2D system.

        Args:
            ssfi (interp1d): Interpolator for the structure factor.
            qcap (float): Cutoff value for the computation.
            number_of_particles (int): Number of particles in the system.

        Returns:
            float: Computed discrete correction value.
        """
        _lambda = 1.0 / np.sqrt(2.0)
        arg_const = np.sqrt(2.0 * np.pi) / np.sqrt(number_of_particles)
        lmax = int(np.ceil(qcap / arg_const / np.sqrt(2.0)))
        pref = 1.0 / (np.sqrt(8.0 * np.pi) * _lambda * np.sqrt(number_of_particles))
        lvals = np.arange(0, lmax + 1, dtype=int)
        L1, L2 = np.meshgrid(lvals, lvals, indexing="ij")
        lnorm = np.sqrt(L1.astype(float) ** 2 + L2.astype(float) ** 2)
        non_zero = (L1 > 0).astype(int) + (L2 > 0).astype(int)
        mult = np.zeros_like(lnorm, dtype=float)
        mult[non_zero == 2] = 4.0  # ++
        mult[non_zero == 1] = 2.0  # +0
        q = arg_const * lnorm
        with np.errstate(divide="ignore", invalid="ignore"):
            term = (ssfi(q) - 1.0) / np.where(lnorm == 0.0, np.inf, lnorm) * mult
            term[np.isinf(term) | np.isnan(term)] = 0.0
        return pref * float(np.sum(term))

    def _compute_discrete_3D(
        self, ssfi: interp1d, qcap: float, number_of_particles: int
    ) -> float:
        """
        Computes the discrete part of the correction for a 3D system.

        Args:
            ssfi (interp1d): Interpolator for the structure factor.
            qcap (float): Cutoff value for the computation.
            number_of_particles (int): Number of particles in the system.

        Returns:
            float: Computed discrete correction value.
        """
        _lambda = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
        arg_const = (8.0 * np.pi / 3.0) ** (1.0 / 3.0) / (
            number_of_particles ** (1.0 / 3.0)
        )
        lmax = int(np.ceil(qcap / arg_const / np.sqrt(3.0)))
        pref = (
            (2.0 / (3.0 * np.pi))
            * ((3.0 / (8.0 * np.pi)) ** (2.0 / 3.0))
            * (1.0 / _lambda)
            * (number_of_particles ** (-1.0 / 3.0))
        )
        total = 0.0
        l2 = np.arange(0, lmax + 1, dtype=int)
        l3 = np.arange(0, lmax + 1, dtype=int)
        L2, L3 = np.meshgrid(l2, l3, indexing="ij")
        for l1 in range(0, lmax + 1):
            L2f, L3f = L2.astype(float), L3.astype(float)
            lnorm = np.sqrt((l1 * l1) + L2f * L2f + L3f * L3f)
            non_zero = (
                (1 if l1 > 0 else 0) + (L2f > 0).astype(int) + (L3f > 0).astype(int)
            )
            mult = np.zeros_like(lnorm, dtype=float)
            mult[non_zero == 3] = 8.0  # +++
            mult[non_zero == 2] = 4.0  # ++0
            mult[non_zero == 1] = 2.0  # +00
            q = arg_const * lnorm
            with np.errstate(divide="ignore", invalid="ignore"):
                term = (ssfi(q) - 1.0) / (lnorm**2) * mult
                term[np.isinf(term) | np.isnan(term)] = 0.0
            if l1 == 0:
                term[0, 0] = 0.0
            total += float(np.sum(term))
        return pref * total

    def _compute_madelung(self, number_of_particles):
        """
        Computes the Madelung constant for the system.

        Args:
            number_of_particles (int): Number of particles in the system.

        Returns:
            float: Computed Madelung constant.
        """
        return (
            -3.90026492 * (1.0 / np.sqrt(np.pi)) * (number_of_particles ** (-1.0 / 2.0))
            if self._is_2D()
            else -2.837297
            * ((3.0 / (4.0 * np.pi)) ** (1.0 / 3.0))
            * (number_of_particles ** (-1.0 / 3.0))
        )
