from __future__ import annotations
import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from . import database
from . import hf
from .dimension import Dimension
from .output import DataBase


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
        self.db_handler = database.DataBaseHandler()
        # Let's postpone specifying the database name to later work when we enable it for
        # all schemes as well
        self.results: Result = Result()
        self.solver: hf.Solver | None = None
        self.inputs: hf.Input | None = None
        self.drs: float | None = None
        self.runs: dict[int, Run] = {}
        self.correction = Correction()
        self.target_run_id: int | None = None

    def compute(
        self,
        solver: hf.Solver,
        inputs: hf.Input,
        drs: float,
        number_of_particles: list[int],
    ) -> int:
        # Assign member variables
        self.solver = solver
        self.inputs = inputs
        self.drs = drs
        self._solve_scheme_for_all_coupling()
        self.correction.compute(self.runs, self.inputs, number_of_particles)
        return self.target_run_id
        # # Write these keys to the scheme's *target* run
        # self.db_handler.run_id = self.target_run_id
        # self.db_handler.insert_results(
        #     {
        #         "fsc_uint": np.asarray(self.correction.internal_energy, dtype=float),
        #         "fsc_fxc": np.asarray(self.correction.free_energy, dtype=float),
        #     },
        #     conflict_mode=database.DataBaseHandler.ConflictMode.UPDATE,
        # )
        # # Mirror into in-memory Result
        # self.results.internal_energy = np.asarray(self.correction.internal_energy, dtype=float)
        # self.results.free_energy = np.asarray(self.correction.free_energy, dtype=float)
        # self.results.run_id = self.target_run_id
        # return self.target_run_id

    def _build_rs_grid(self) -> np.ndarray:
        if self.drs <= 0:
            raise ValueError(
                f"Invalid coupling resolution, it must be at greater than 0"
            )
        target = self.inputs.coupling
        grid = np.arange(0.0, target + self.drs, self.drs)
        if len(grid) == 0 or not math.isclose(
            grid[-1], target, rel_tol=0, abs_tol=FLOAT_TOLERANCE
        ):
            grid = np.append(grid, target)
        if not math.isclose(grid[0], 0.0, rel_tol=0, abs_tol=FLOAT_TOLERANCE):
            grid = np.insert(grid, 0, 0.0)
        return grid.astype(float)

    def _solve_scheme_for_all_coupling(self):
        rs_grid = self._build_rs_grid()
        self._search_results_in_the_database(rs_grid)
        filtered_rs_grid = [
            rs for rs in rs_grid if self._key_from_float(rs) not in self.runs
        ]
        target_coupling = self.inputs.coupling
        for rs in filtered_rs_grid:
            print(f"--- Computing run for rs={rs}")
            self.inputs.coupling = rs
            self.solver.compute(self.inputs)
            key = self._key_from_float(rs)
            self.runs[key] = Run(rs, self.solver.run_id)
        self.inputs.coupling = target_coupling
        target_key = self._key_from_float(target_coupling)
        self.target_run_id = self.runs[target_key].id
        print(f"--- Runs completed. Target run ID: {self.target_run_id}")

    def _key_from_float(self, x):
        precision = int(np.abs(np.log10(self.drs)))
        return int(round(x * 10**precision))

    def _search_results_in_the_database(self, rs_grid):
        runs = self.db_handler.inspect_runs()
        database_keys = database.DataBaseHandler.TableKeys
        filtered_runs = [
            run
            for run in runs
            if run[database_keys.THEORY.value] == self.inputs.theory
            and run[database_keys.DEGENERACY.value] == self.inputs.degeneracy
        ]
        for rs in rs_grid:
            # It's better not to use floats as keys in a dictionary
            key = self._key_from_float(rs)
            run = next(
                (
                    run
                    for run in filtered_runs
                    if self._key_from_float(run[database_keys.COUPLING.value]) == key
                ),
                None,
            )
            if run is None:
                continue
            run_id = run[database_keys.PRIMARY_KEY.value]
            run_inputs = self.db_handler.get_inputs(run_id)
            if (
                run_inputs["cutoff"] == self.inputs.cutoff
                and run_inputs["matsubara"] == self.inputs.matsubara
                and run_inputs["resolution"] == self.inputs.resolution
            ):
                self.runs[key] = Run(rs, run_id)


@dataclass
class Result:
    internal_energy: np.ndarray | None = None
    free_energy: np.ndarray | None = None
    run_id: int | None = None


class Correction:

    def __init__(self):
        self.continuous = ContinuousCorrection()
        self.internal_energy: float | None = None
        self.free_energy: float | None = None

    def compute(
        self,
        runs: dict[float, Run],
        inputs: hf.Input,
        number_of_particles: int,
    ):
        self.continuous.compute(runs, inputs)
        interp_kind = "cubic" if len(runs) >= 4 else "linear"
        rs_grid = []
        for _, run in runs.items():
            rs_grid.append(run.rs)
        rs_grid.sort()
        rs_uint_vals = []
        target_rs_uint: float | None = None
        for key, run in runs.items():
            continuous = self.continuous.results[key]
            ssfi = interp1d(
                self.continuous.wvg,
                self.continuous.ssf,
                kind="cubic",
                fill_value="extrapolate",
            )
            discrete = self._compute_result(ssfi, inputs.cutoff, number_of_particles)
            madelung = self._compute_madelung(number_of_particles)
            rs_uint = continuous - discrete - 0.5 * madelung
            rs_uint_vals.append(rs_uint)
            if math.isclose(
                run.rs,
                inputs.coupling,
                rel_tol=0,
                abs_tol=FLOAT_TOLERANCE,
            ):
                target_rs_uint = rs_uint
        if target_rs_uint is None:
            raise RuntimeError("Target rs not found in grid")
        # uint FSC at target
        self.internal_energy = target_rs_uint / inputs.coupling
        # fxc FSC: integrate rs*uint(rs) from 0 to target_rs, then divide by rs^2
        rs_grid_np = np.asarray(rs_grid, dtype=float)
        rs_uint_np = np.asarray(rs_uint_vals, dtype=float)
        interp = interp1d(
            rs_grid_np, rs_uint_np, kind=interp_kind, fill_value="extrapolate"
        )
        integral_val = quad(
            lambda r: float(interp(r)), 0.0, float(inputs.coupling), limit=200
        )[0]
        self.free_energy = integral_val / (inputs.coupling**2.0)

    def _compute_result(
        self, ssfi: interp1d, cutoff: float, number_of_particles: int
    ) -> float:
        return (
            self._compute_result_2D(ssfi, cutoff, number_of_particles)
            if self.continuous.is_2D
            else self._compute_result_3D(ssfi, cutoff, number_of_particles)
        )

    def _compute_result_2D(
        ssfi: interp1d, qcap: float, number_of_particles: int
    ) -> float:
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

    def _compute_result_3D(
        self, ssfi: interp1d, qcap: float, number_of_particles: int
    ) -> float:
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
        return (
            -3.90026492 * (1.0 / np.sqrt(np.pi)) * (number_of_particles ** (-1.0 / 2.0))
            if self.continuous.is_2D
            else -2.837297
            * ((3.0 / (4.0 * np.pi)) ** (1.0 / 3.0))
            * (number_of_particles ** (-1.0 / 3.0))
        )


class ContinuousCorrection:

    def __init__(self):
        self.results: dict[int, float] = {}
        self.is_2D = False
        self.ssf: None | list[float] = None
        self.wvg: None | list[float] = None

    def compute(self, runs: dict[float, Run], inputs: hf.Input):
        if not self.results:
            self.is_2D = inputs.dimension == Dimension._2D
            for key, run in runs.items():
                self._read_wvg_ssf(run.id)
                ssfi = interp1d(
                    self.wvg, self.ssf, kind="cubic", fill_value="extrapolate"
                )
                self.results[key] = self._compute_result(ssfi, cutoff=inputs.cutoff)

    def _compute_result(self, ssfi: interp1d, cutoff: float) -> float:
        return (
            self._compute_result_2D(ssfi, cutoff)
            if self.is_2D
            else self._compute_result_3D(ssfi, cutoff)
        )

    def _compute_result_3D(self, ssfi: interp1d, cutoff: float) -> float:
        _lambda = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
        val = quad(lambda q: float(ssfi(q) - 1.0), 0.0, cutoff, limit=200)[0]
        return val / (np.pi * _lambda)

    def _compute_result_2D(self, ssfi: interp1d, cutoff: float) -> float:
        _lambda = 1.0 / np.sqrt(2.0)
        val = quad(lambda q: float(ssfi(q) - 1.0), 0.0, cutoff, limit=200)[0]
        return val / (2.0 * _lambda)

    def _read_wvg_ssf(self, run_id: int) -> None:
        data = DataBase.read_run(int(run_id))
        res = data.get("results", {})
        self.wvg = res.get("wvg", None)
        self.ssf = res.get("ssf", None)
        if self.wvg is None or self.ssf is None:
            raise ValueError(f"Malformed results for Run {run_id}.")
