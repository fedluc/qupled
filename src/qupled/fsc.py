from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from . import database
from . import hf
from .dimension import Dimension
from .output import DataBase


class Solver:
    """
    Post-processing FSC solver. It computes and writes results
    into the *scheme* run at (rs=target_rs, theta, dimension).
    """

    def __init__(
        self,
        *,
        scheme_module,
        scheme_kwargs: Optional[dict] = None,
        database_name: str | None = None,
    ):
        self.scheme_module = scheme_module
        self.scheme_kwargs = dict(scheme_kwargs or {})
        self.dbh = database.DataBaseHandler(database_name)
        self.results: Result = Result()

    def compute(self, inputs: "Input", N_values: Iterable[int]) -> int:
        # Dimension handling via enum (default to 3D)
        dim_enum: Dimension = (
            inputs.dimension
            if isinstance(inputs.dimension, Dimension)
            else Dimension._3D
        )
        is_2d = dim_enum is Dimension._2D

        # If user passed a dimension in scheme_kwargs, ensure consistency
        kw_dim = self.scheme_kwargs.get("dimension", None)
        if kw_dim is not None:
            if not isinstance(kw_dim, Dimension):
                raise TypeError(
                    "scheme_kwargs['dimension'] must be a qupled.dimension.Dimension enum."
                )
            if kw_dim is not dim_enum:
                raise ValueError(
                    f"Dimension mismatch: inputs.dimension={dim_enum} vs scheme_kwargs['dimension']={kw_dim}"
                )

        scheme_kwargs_dim = dict(self.scheme_kwargs)
        scheme_kwargs_dim["dimension"] = dim_enum

        theory_tag = _infer_theory_tag(
            self.scheme_module,
            rs=inputs.coupling,
            theta=inputs.degeneracy,
            **scheme_kwargs_dim,
        )

        # Build rs-grid
        rs_grid = _build_rs_grid(inputs.coupling, inputs.drs)

        # Ensure runs for every rs on the grid
        run_ids_by_rs: Dict[float, int] = {}
        for rs in rs_grid:
            rid = _ensure_run_for(
                self.dbh,
                self.scheme_module,
                theory_tag,
                float(rs),
                float(inputs.degeneracy),
                dim_enum,
                **scheme_kwargs_dim,
            )
            run_ids_by_rs[float(rs)] = rid

        target_run_id = run_ids_by_rs[float(rs_grid[-1])]

        # Set λ per dimension
        lambda_val = (
            (1.0 / np.sqrt(2.0)) if is_2d else (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
        )

        # Interp kind
        interp_kind = "cubic" if len(rs_grid) >= 4 else "linear"

        # Precompute continuous term at each rs (N-independent)
        cont_by_rs: Dict[float, float] = {}
        qcap_by_rs: Dict[float, float] = {}
        for rs in rs_grid:
            wvg, ssf = _read_wvg_ssf(run_ids_by_rs[float(rs)])
            S_interp = _make_S_interp(wvg, ssf)
            qcap = min(inputs.cutoff, float(wvg[-1]))
            qcap_by_rs[float(rs)] = qcap
            cont_by_rs[float(rs)] = (
                _continuous_term_2d(S_interp, qcap, lambda_val)
                if is_2d
                else _continuous_term_3d(S_interp, qcap, lambda_val)
            )

        fsc_uint_list: List[float] = []
        fsc_fxc_list: List[float] = []

        for N in sorted(N_values):
            rs_uint_vals: List[float] = []
            target_rs_uint: Optional[float] = None

            for rs in rs_grid:
                wvg, ssf = _read_wvg_ssf(run_ids_by_rs[float(rs)])
                S_interp = _make_S_interp(wvg, ssf)
                qcap = qcap_by_rs[float(rs)]
                continuous = cont_by_rs[float(rs)]

                if is_2d:
                    discrete = _discrete_term_2d(N, qcap, S_interp, lambda_val)
                    madelung = _madelung_2d(N)
                else:
                    discrete = _discrete_term_3d(N, qcap, S_interp, lambda_val)
                    madelung = _madelung_3d(N)

                rs_uint = continuous - discrete - 0.5 * madelung
                rs_uint_vals.append(rs_uint)

                if math.isclose(
                    float(rs), float(inputs.coupling), rel_tol=0, abs_tol=1e-12
                ):
                    target_rs_uint = rs_uint

            if target_rs_uint is None:
                raise RuntimeError("Target rs not found in grid")

            # uint FSC at target
            fsc_uint = target_rs_uint / float(inputs.coupling)
            fsc_uint_list.append(float(fsc_uint))

            # fxc FSC: integrate rs*uint(rs) from 0 to target_rs, then divide by rs^2
            rs_grid_np = np.asarray(rs_grid, dtype=float)
            rs_uint_np = np.asarray(rs_uint_vals, dtype=float)
            interp = interp1d(
                rs_grid_np, rs_uint_np, kind=interp_kind, fill_value="extrapolate"
            )
            integral_val = quad(
                lambda r: float(interp(r)), 0.0, float(inputs.coupling), limit=200
            )[0]
            fsc_fxc = integral_val / (float(inputs.coupling) ** 2.0)
            fsc_fxc_list.append(float(fsc_fxc))

        # Write these keys to the scheme's *target* run
        self.dbh.run_id = target_run_id
        self.dbh.insert_results(
            {
                "fsc_uint": np.asarray(fsc_uint_list, dtype=float),
                "fsc_fxc": np.asarray(fsc_fxc_list, dtype=float),
            },
            conflict_mode=database.DataBaseHandler.ConflictMode.UPDATE,
        )

        # Mirror into in-memory Result
        self.results.fsc_uint = np.asarray(fsc_uint_list, dtype=float)
        self.results.fsc_fxc = np.asarray(fsc_fxc_list, dtype=float)
        self.results.run_id = target_run_id
        return target_run_id


@dataclass
class Input(hf.Input):
    drs: float = 0.1


@dataclass
class Result:
    fsc_uint: np.ndarray | None = None
    fsc_fxc: np.ndarray | None = None
    run_id: int | None = None


# ============================ Internal helpers =============================


def _build_rs_grid(target_rs: float, drs: float) -> np.ndarray:
    if drs <= 0:
        raise ValueError("drs must be > 0")
    grid = np.arange(0.0, target_rs + drs, drs)
    if len(grid) == 0 or not math.isclose(
        grid[-1], target_rs, rel_tol=0, abs_tol=1e-12
    ):
        grid = np.append(grid, target_rs)
    if not math.isclose(grid[0], 0.0, rel_tol=0, abs_tol=1e-15):
        grid = np.insert(grid, 0, 0.0)
    return grid.astype(float)


def _infer_theory_tag(scheme_module, rs: float, theta: float, **scheme_kwargs) -> str:
    tmp = scheme_module.Input(rs, theta, **scheme_kwargs)
    return str(getattr(tmp, "theory", scheme_module.__name__.split(".")[-1].upper()))


def _dim_value(dim_enum: Dimension) -> str:
    """Return '2D' or '3D' from a Dimension enum."""
    return dim_enum.value


def _ensure_run_for(
    dbh: database.DataBaseHandler,
    scheme_module,
    theory: str,
    rs: float,
    theta: float,
    dim_enum: Dimension,
    **scheme_kwargs,
) -> int:
    rid = _find_run_id(dbh, theory, rs, theta, dim_enum)
    if rid is not None:
        return rid

    kw_no_dim = {k: v for k, v in scheme_kwargs.items() if k != "dimension"}
    inputs = scheme_module.Input(rs, theta, **kw_no_dim)
    inputs.dimension = dim_enum
    solver = scheme_module.Solver()
    solver.compute(inputs)
    if getattr(solver, "run_id", None) is None:
        raise RuntimeError("Scheme solver did not expose run_id after compute().")

    run_data = DataBase.read_run(int(solver.run_id))
    run_dim_val = Dimension.from_dict(
        run_data["inputs"]["dimension"]
    ).value  # "2D"/"3D"
    if run_dim_val != _dim_value(dim_enum):
        raise RuntimeError(
            f"Computed run has dimension {run_dim_val}, expected {_dim_value(dim_enum)}."
        )
    return int(solver.run_id)


def _find_run_id(
    dbh: database.DataBaseHandler,
    theory: str,
    rs: float,
    theta: float,
    dim_enum: Dimension,
) -> Optional[int]:
    """
    Find a run that matches theory, (rs, theta) and desired dimension.
    Uses DataBaseHandler.inspect_runs() to list, then DataBase.read_run() to check inputs.
    """
    desired_val = _dim_value(dim_enum)  # "2D" or "3D"
    runs = dbh.inspect_runs()
    candidates = []
    for r in runs:
        if r.get("theory") != theory:
            continue
        if "coupling" not in r or "degeneracy" not in r:
            continue
        if not (
            math.isclose(float(r["coupling"]), rs, rel_tol=0, abs_tol=1e-12)
            and math.isclose(float(r["degeneracy"]), theta, rel_tol=0, abs_tol=1e-12)
        ):
            continue
        run_full = DataBase.read_run(int(r["id"]))
        dim_val = Dimension.from_dict(
            run_full["inputs"]["dimension"]
        ).value  # "2D"/"3D"
        if dim_val == desired_val:
            candidates.append(r)
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x.get("date", ""), x.get("time", "")))
    return int(candidates[-1]["id"])


def _read_wvg_ssf(run_id: int) -> Tuple[np.ndarray, np.ndarray]:
    data = DataBase.read_run(int(run_id))
    res = data.get("results", {})
    wvg = res.get("wvg", None)
    ssf = res.get("ssf", None)
    if wvg is None or ssf is None:
        raise ValueError(f"Run {run_id} lacks 'wvg'/'ssf' results required for FSC.")
    return np.asarray(wvg, dtype=float), np.asarray(ssf, dtype=float)


def _make_S_interp(wvg: np.ndarray, ssf: np.ndarray):
    return interp1d(
        wvg,
        ssf,
        kind="cubic",
        bounds_error=False,
        assume_sorted=False,
    )


# =========================== Continuous terms ============================


def _continuous_term_3d(S_interp, qcap, lambda_val):
    val = quad(lambda q: float(S_interp(q) - 1.0), 0.0, float(qcap), limit=200)[0]
    return val / (np.pi * lambda_val)


def _continuous_term_2d(S_interp, qcap, lambda_val):
    val = quad(lambda q: float(S_interp(q) - 1.0), 0.0, float(qcap), limit=200)[0]
    return val / (2.0 * lambda_val)


# ============================ Discrete terms =============================


def _discrete_term_3d(N: int, qcap: float, S_interp, lambda_val: float) -> float:
    """
    Octant sum with multiplicities: +++, ++0, +00.
    """
    arg_const = (8.0 * np.pi / 3.0) ** (1.0 / 3.0) / (N ** (1.0 / 3.0))
    lmax = int(np.ceil(qcap / arg_const / np.sqrt(3.0)))
    pref = (
        (2.0 / (3.0 * np.pi))
        * ((3.0 / (8.0 * np.pi)) ** (2.0 / 3.0))
        * (1.0 / lambda_val)
        * (N ** (-1.0 / 3.0))
    )

    total = 0.0
    l2 = np.arange(0, lmax + 1, dtype=int)
    l3 = np.arange(0, lmax + 1, dtype=int)
    L2, L3 = np.meshgrid(l2, l3, indexing="ij")

    for l1 in range(0, lmax + 1):
        L2f, L3f = L2.astype(float), L3.astype(float)
        lnorm = np.sqrt((l1 * l1) + L2f * L2f + L3f * L3f)

        non_zero = (1 if l1 > 0 else 0) + (L2f > 0).astype(int) + (L3f > 0).astype(int)
        mult = np.zeros_like(lnorm, dtype=float)
        mult[non_zero == 3] = 8.0  # +++
        mult[non_zero == 2] = 4.0  # ++0
        mult[non_zero == 1] = 2.0  # +00

        q = arg_const * lnorm
        with np.errstate(divide="ignore", invalid="ignore"):
            term = (S_interp(q) - 1.0) / (lnorm**2) * mult
            term[np.isinf(term) | np.isnan(term)] = 0.0
        if l1 == 0:
            term[0, 0] = 0.0
        total += float(np.sum(term))
    return pref * total


def _discrete_term_2d(N: int, qcap: float, S_interp, lambda_val: float) -> float:
    arg_const = np.sqrt(2.0 * np.pi) / np.sqrt(N)
    lmax = int(np.ceil(qcap / arg_const / np.sqrt(2.0)))
    pref = 1.0 / (np.sqrt(8.0 * np.pi) * lambda_val * np.sqrt(N))

    total = 0.0
    lvals = np.arange(0, lmax + 1, dtype=int)
    L1, L2 = np.meshgrid(lvals, lvals, indexing="ij")
    lnorm = np.sqrt(L1.astype(float) ** 2 + L2.astype(float) ** 2)

    non_zero = (L1 > 0).astype(int) + (L2 > 0).astype(int)
    mult = np.zeros_like(lnorm, dtype=float)
    mult[non_zero == 2] = 4.0  # ++
    mult[non_zero == 1] = 2.0  # +0

    q = arg_const * lnorm
    with np.errstate(divide="ignore", invalid="ignore"):
        term = (S_interp(q) - 1.0) / np.where(lnorm == 0.0, np.inf, lnorm) * mult
        term[np.isinf(term) | np.isnan(term)] = 0.0
    return pref * float(np.sum(term))


# ============================ Madelung terms =============================


def _madelung_3d(N: int) -> float:
    return -2.837297 * ((3.0 / (4.0 * np.pi)) ** (1.0 / 3.0)) * (N ** (-1.0 / 3.0))


def _madelung_2d(N: int) -> float:
    return -3.90026492 * (1.0 / np.sqrt(np.pi)) * (N ** (-1.0 / 2.0))
