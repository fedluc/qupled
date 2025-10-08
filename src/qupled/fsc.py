from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from . import database
from . import hf 

class Solver:
    """
    Post-processing FSC solver. It computes and writes results 
    into the *scheme* run at (rs=target_rs, theta, dimension).
    """

    def __init__(self, *, scheme_module, scheme_kwargs: Optional[dict] = None, database_name: str | None = None):
        """
        Args:
            scheme_module: module providing Solver and Input (e.g., qupled.stls)
            scheme_kwargs: kwargs passed to the scheme's Input/Solver (must include 'dimension': 'D2'|'D3')
            database_name: optional DB name
        """
        self.scheme_module = scheme_module
        self.scheme_kwargs = scheme_kwargs or {}
        self.dbh = database.DataBaseHandler(database_name)
        self.results: Result = Result()

    def compute(self, inputs: "Input", N_values: Iterable[int]) -> int:
        dim = (str(inputs.dimension).upper() or "D3")
        if dim not in {"D2", "D3"}:
            raise ValueError(f"Unsupported dimension: {inputs.dimension!r}")

        kw_dim = str(self.scheme_kwargs.get("dimension", "")).upper()
        if kw_dim and kw_dim != dim:
            raise ValueError(
                f"Dimension mismatch: inputs.dimension={inputs.dimension} vs "
                f"scheme_kwargs['dimension']={kw_dim}"
            )

        scheme_kwargs_dim = dict(self.scheme_kwargs)
        scheme_kwargs_dim.setdefault("dimension", dim)

        theory_tag = _infer_theory_tag(
            self.scheme_module, rs=inputs.coupling, theta=inputs.degeneracy, **scheme_kwargs_dim
        )

        rs_grid = _build_rs_grid(inputs.coupling, inputs.drs)

        # ensure runs (dimension-matching) for each rs in the grid
        run_ids_by_rs: Dict[float, int] = {}
        for rs in rs_grid:
            rid = _ensure_run_for(
                self.dbh, self.scheme_module, theory_tag, float(rs), float(inputs.degeneracy),
                desired_dim=dim, **scheme_kwargs_dim
            )
            run_ids_by_rs[float(rs)] = rid

        target_run_id = run_ids_by_rs[float(rs_grid[-1])]

        # set lambda per dimension
        if dim == "D2":
            lambda_val = 1.0 / np.sqrt(2.0)
        else: 
            lambda_val = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)

        interp_kind = "cubic" if len(rs_grid) >= 4 else "linear"

        # precompute continuous term at each rs (N-independent)
        cont_by_rs: Dict[float, float] = {}
        qcap_by_rs: Dict[float, float] = {}
        for rs in rs_grid:
            wvg, ssf = _read_wvg_ssf(self.dbh, run_ids_by_rs[float(rs)])
            S_interp = _make_S_interp(wvg, ssf)
            qcap = min(inputs.cutoff, float(wvg[-1]))
            qcap_by_rs[float(rs)] = qcap
            cont_by_rs[float(rs)] = (
                _continuous_term_2d(S_interp, qcap, lambda_val)
                if dim == "D2" else
                _continuous_term_3d(S_interp, qcap, lambda_val)
            )

        fsc_uint_list: List[float] = []
        fsc_fxc_list:  List[float] = []

        for N in sorted(N_values):
            rs_uint_vals: List[float] = []
            target_rs_uint: Optional[float] = None

            for rs in rs_grid:
                wvg, ssf = _read_wvg_ssf(self.dbh, run_ids_by_rs[float(rs)])
                S_interp = _make_S_interp(wvg, ssf)
                qcap = qcap_by_rs[float(rs)]
                continuous = cont_by_rs[float(rs)]

                if dim == "D2":
                    discrete = _discrete_term_2d(N, qcap, S_interp, lambda_val)
                    madelung = _madelung_2d(N)
                else:  # "D3"
                    discrete = _discrete_term_3d(N, qcap, S_interp, lambda_val)
                    madelung = _madelung_3d(N)

                rs_uint = continuous - discrete - 0.5 * madelung 
                rs_uint_vals.append(rs_uint)

                if math.isclose(float(rs), float(inputs.coupling), rel_tol=0, abs_tol=1e-12):
                    target_rs_uint = rs_uint

            if target_rs_uint is None:
                raise RuntimeError("Target rs not found in grid")

            fsc_uint = target_rs_uint / float(inputs.coupling)
            fsc_uint_list.append(float(fsc_uint))

            rs_grid_np = np.asarray(rs_grid, dtype=float)
            rs_uint_np = np.asarray(rs_uint_vals, dtype=float)
            interp = interp1d(rs_grid_np, rs_uint_np, kind=interp_kind, fill_value="extrapolate")
            integral_val = quad(lambda r: float(interp(r)), 0.0, float(inputs.coupling), limit=200)[0]
            fsc_fxc = integral_val / (float(inputs.coupling) ** 2.0)
            fsc_fxc_list.append(float(fsc_fxc))

        self.dbh.run_id = target_run_id
        self.dbh.insert_results(
            {"fsc_uint": np.asarray(fsc_uint_list, dtype=float),
            "fsc_fxc":  np.asarray(fsc_fxc_list,  dtype=float)},
            conflict_mode=database.DataBaseHandler.ConflictMode.UPDATE,
        )

        self.results.fsc_uint = np.asarray(fsc_uint_list, dtype=float)
        self.results.fsc_fxc  = np.asarray(fsc_fxc_list,  dtype=float)
        self.results.run_id   = target_run_id
        return target_run_id


@dataclass
class Input(hf.Input):
    """
    FSC post-processing inputs. Inherit all standard fields from hf.Input:
    """
    drs: float = 0.1


@dataclass
class Result:
    fsc_uint: np.ndarray | None = None
    fsc_fxc:  np.ndarray | None = None
    run_id:   int | None = None


# --------------------------- Internal helpers --------------------------------

def _build_rs_grid(target_rs: float, drs: float) -> np.ndarray:
    if drs <= 0:
        raise ValueError("drs must be > 0")
    grid = np.arange(0.0, target_rs + drs, drs)
    if len(grid) == 0 or not math.isclose(grid[-1], target_rs, rel_tol=0, abs_tol=1e-12):
        grid = np.append(grid, target_rs)
    if not math.isclose(grid[0], 0.0, rel_tol=0, abs_tol=1e-15):
        grid = np.insert(grid, 0, 0.0)
    return grid.astype(float)

def _infer_theory_tag(scheme_module, rs: float, theta: float, **scheme_kwargs) -> str:
    tmp = scheme_module.Input(rs, theta, **scheme_kwargs)
    return str(getattr(tmp, "theory", scheme_module.__name__.split(".")[-1].upper()))

def _ensure_run_for(dbh: database.DataBaseHandler, scheme_module, theory: str,
                    rs: float, theta: float, desired_dim: str, **scheme_kwargs) -> int:
    rid = _find_run_id(dbh, theory, rs, theta, desired_dim)
    if rid is not None:
        return rid
    # Compute run with the provided kwargs
    solver = scheme_module.Solver()
    inputs = scheme_module.Input(rs, theta, **scheme_kwargs)
    solver.compute(inputs)
    if getattr(solver, "run_id", None) is None:
        raise RuntimeError("Scheme solver did not expose run_id after compute().")
    # confirm dimension matches
    inp = dbh.get_inputs(int(solver.run_id), names=["dimension"])
    new_dim = str(inp.get("dimension", "")).upper()
    if new_dim != desired_dim:
        raise RuntimeError(f"Computed run has dimension {new_dim}, expected {desired_dim}.")
    return int(solver.run_id)

def _find_run_id(dbh: database.DataBaseHandler, theory: str, rs: float, theta: float,
                 desired_dim: str) -> Optional[int]:
    runs = dbh.inspect_runs()
    candidates = []
    for r in runs:
        if r.get("theory") != theory:
            continue
        if "coupling" not in r or "degeneracy" not in r:
            continue
        if not (math.isclose(float(r["coupling"]), rs, rel_tol=0, abs_tol=1e-12)
                and math.isclose(float(r["degeneracy"]), theta, rel_tol=0, abs_tol=1e-12)):
            continue
        inp = dbh.get_inputs(int(r["id"]), names=["dimension"])
        run_dim = str(inp.get("dimension", "D3")).upper()  # default to D3 if missing
        if run_dim == desired_dim:
            candidates.append(r)
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x.get("date", ""), x.get("time", "")))
    return int(candidates[-1]["id"])

def _read_wvg_ssf(dbh: database.DataBaseHandler, run_id: int) -> Tuple[np.ndarray, np.ndarray]:
    out = dbh.get_results(run_id, names=["wvg", "ssf"])
    wvg = out.get("wvg", None)
    ssf = out.get("ssf", None)
    if wvg is None or ssf is None:
        raise ValueError(f"Run {run_id} lacks 'wvg'/'ssf' results required for FSC.")
    return np.asarray(wvg, dtype=float), np.asarray(ssf, dtype=float)

def _make_S_interp(wvg: np.ndarray, ssf: np.ndarray):
    return interp1d(
        wvg, ssf, kind="cubic",
        bounds_error=False,
        #fill_value=(1.0, ssf[-1]),  # solver should cap q at cutoff; tail shouldn’t be used
        assume_sorted=False,
    )

# --------------------------- Continuous terms --------------------------------

def _continuous_term_3d(S_interp, qcap, lambda_val):
    val = quad(lambda q: float(S_interp(q) - 1.0), 0.0, float(qcap), limit=200)[0]
    return val / (np.pi * lambda_val)

def _continuous_term_2d(S_interp, qcap, lambda_val):
    val = quad(lambda q: float(S_interp(q) - 1.0), 0.0, float(qcap), limit=200)[0]
    return val / (2.0 * lambda_val)

# --------------------------- Discrete terms --------------------------------

def _discrete_term_3d(N: int, qcap: float, S_interp, lambda_val: float) -> float:
    """
    The sum over the discrete reciprocal vectors is computed only for 
    positive values of l1,l2,l3 due to the negative values being symmetric, i.e.
    we care only about the +++ octant. Then checks are made for how many zeros 
    are present in each combination. In the case that only one of 
    the vectors is non-zero (+00) then the resulting value is multiplied by 2, 
    if two are non-zero (++0) then it's multiplied by 4, and if all three are 
    non-zero (+++) then it's multiplied by 8. Additionally, for further speedup
    we take advantage of NumPy vectorized operations by precomputing the norms 
    for each value of l1, i.e. each 2D slice of the l2,l3 values. 
    """
    arg_const = (8.0 * np.pi / 3.0) ** (1.0/3.0) / (N ** (1.0/3.0))
    lmax = int(np.ceil(qcap / arg_const / np.sqrt(3.0)))
    pref = (2.0 / (3.0 * np.pi)) * ((3.0 / (8.0 * np.pi)) ** (2.0/3.0)) * (1.0 / lambda_val) * (N ** (-1.0/3.0))

    total = 0.0
    l2 = np.arange(0, lmax + 1, dtype=int)
    l3 = np.arange(0, lmax + 1, dtype=int)
    L2, L3 = np.meshgrid(l2, l3, indexing="ij")

    for l1 in range(0, lmax + 1):
        L2f, L3f = L2.astype(float), L3.astype(float)
        lnorm = np.sqrt((l1*l1) + L2f*L2f + L3f*L3f)

        non_zero = (1 if l1 > 0 else 0) + (L2f > 0).astype(int) + (L3f > 0).astype(int)
        mult = np.zeros_like(lnorm, dtype=float)
        mult[non_zero == 3] = 8.0 #+++
        mult[non_zero == 2] = 4.0 #++0
        mult[non_zero == 1] = 2.0 #+00

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
    lnorm = np.sqrt(L1.astype(float)**2 + L2.astype(float)**2)

    non_zero = (L1 > 0).astype(int) + (L2 > 0).astype(int)
    mult = np.zeros_like(lnorm, dtype=float)
    mult[non_zero == 2] = 4.0   # ++
    mult[non_zero == 1] = 2.0   # +0

    q = arg_const * lnorm
    with np.errstate(divide="ignore", invalid="ignore"):
        term = (S_interp(q) - 1.0) / np.where(lnorm == 0.0, np.inf, lnorm) * mult
        term[np.isinf(term) | np.isnan(term)] = 0.0
    return pref * float(np.sum(term))

# --------------------------- Madelung terms --------------------------------

def _madelung_3d(N: int) -> float:
    return -2.837297 * ((3.0 / (4.0 * np.pi)) ** (1.0/3.0)) * (N ** (-1.0/3.0))

def _madelung_2d(N: int) -> float:
    return -3.90026492 * (1.0 / np.sqrt(np.pi)) * (N ** (-1.0/2.0))
