# qupled/itcf.py
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy import integrate, optimize

from . import database
from . import hf
from . import native as qnative


# ============================== Public API ===============================

@dataclass
class Input(hf.Input):
    """
    ITCF post-processing inputs. Inherit all fields from hf.Input:
      - coupling (rs), degeneracy (theta), cutoff, dimension, matsubara, etc.
    New fields:
      - dt: time step for t-grid in [0, 0.5].
    """
    dt: float = 0.01


@dataclass
class Result:
    itcf:   np.ndarray | None = None   # (nx, nt)
    x_grid: np.ndarray | None = None   # (nx,)
    t_grid: np.ndarray | None = None   # (nt,)
    run_id: int | None = None


class Solver:
    """
    Post-processing ITCF solver. It computes F(x,t) and writes it into the
    *scheme* run at (rs, theta, dimension). No new run is created.
    """

    def __init__(self, *, scheme_module, scheme_kwargs: Optional[dict] = None,
                 database_name: str | None = None):
        self.scheme_module   = scheme_module
        self.scheme_kwargs   = dict(scheme_kwargs or {})
        self.dbh             = database.DataBaseHandler(database_name)
        self.results: Result = Result()

    def compute(self, inputs: "Input") -> int:
        # ----- enforce dimension (default D3) -----
        dim = (str(inputs.dimension).upper() or "D3")
        if dim not in {"D2", "D3"}:
            raise ValueError(f"Unsupported dimension: {inputs.dimension!r}")
        kw_dim = str(self.scheme_kwargs.get("dimension", "")).upper()
        if kw_dim and kw_dim != dim:
            raise ValueError(f"Dimension mismatch: inputs.dimension={inputs.dimension} "
                             f"vs scheme_kwargs['dimension']={kw_dim}")

        scheme_kwargs_dim = dict(self.scheme_kwargs)
        scheme_kwargs_dim.setdefault("dimension", dim)

        # ----- locate/ensure the scheme run -----
        theory = _infer_theory_tag(self.scheme_module, rs=inputs.coupling,
                                   theta=inputs.degeneracy, **scheme_kwargs_dim)
        run_id = _ensure_run_for(
            self.dbh, self.scheme_module, theory,
            float(inputs.coupling), float(inputs.degeneracy),
            desired_dim=dim, **scheme_kwargs_dim
        )

        # ----- read scheme outputs -----
        wvg, idr, lfc = _read_wvg_idr_lfc(self.dbh, run_id)  # wvg:(nx,), idr:(nx,nl)|(nl,nx), lfc: flexible
        nx = int(wvg.shape[0])

        # time grid t ∈ [0, 0.5]
        if inputs.dt <= 0:
            raise ValueError("dt must be > 0")
        t_grid = np.arange(0.0, 0.5 + inputs.dt, inputs.dt, dtype=float)
        t_grid[-1] = min(t_grid[-1], 0.5)  # clip numeric drift

        # normalize shapes: idr->(nx,nl); lfc->(nx,nl)
        idr, nl = _normalize_idr(idr, nx)
        lfc_eff = _normalize_lfc(lfc, nx, nl)

        theta = float(inputs.degeneracy)
        rs    = float(inputs.coupling)

        # ----- dimension-specific pieces -----
        if dim == "D2":
            # norm = Θ, ip(x) = √2 rs / x
            FHF_mat = _F_HF_2d_matrix(wvg, t_grid, theta)  # (nx, nt)
            norm = theta
            ip_x = np.zeros_like(wvg, dtype=float)
            mask = wvg > 0.0
            ip_x[mask] = np.sqrt(2.0) * rs / wvg[mask]
        else:
            # norm = 3/2 Θ, ip(x) = 4 λ rs / (π x^2), λ = (4/9π)^{1/3}
            lam3d = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
            mu3d  = _mu_3d(theta)
            FHF_mat = _F_HF_3d_matrix(wvg, t_grid, theta, mu3d)
            norm = 1.5 * theta
            ip_x = np.zeros_like(wvg, dtype=float)
            mask = wvg > 0.0
            ip_x[mask] = (4.0 * lam3d * rs) / (np.pi * (wvg[mask] ** 2))

        # ----- Matsubara weights and cosines, including l=0 -----
        l_idx = np.arange(nl, dtype=float)              # 0..nl-1
        w_l   = np.ones(nl, dtype=float)               # w_0=1, w_l=2 for l>=1
        if nl >= 1:
            w_l[1:] = 2.0
        cos_lt = np.cos(2.0 * np.pi * np.outer(l_idx, t_grid))  # (nl, nt)

        # ----- series coefficients a(x,l) -----
        one_m_lfc = 1.0 - lfc_eff                         # (nx, nl)
        denom = 1.0 + (ip_x[:, None] * one_m_lfc * idr)   # (nx, nl)
        denom = np.where(np.abs(denom) < 1e-14, 1e-14, denom)
        a_xl  = (idr * idr) * one_m_lfc / denom           # (nx, nl)

        # weighted sum over Matsubara l
        S_xt = (a_xl * w_l[None, :]) @ cos_lt             # (nx, nt)

        # F(x,t) = F_HF(x,t) - norm * ip(x) * S(x,t)
        F_xt = FHF_mat - (norm * ip_x[:, None]) * S_xt

        # enforce F(0,t) = 0
        F_xt[wvg == 0.0, :] = 0.0

        # ----- write back to DB -----
        self.dbh.run_id = run_id
        self.dbh.insert_results(
            {"itcf": F_xt.astype(float),
             "x_grid": wvg.astype(float),
             "t_grid": t_grid.astype(float)},
            conflict_mode=database.DataBaseHandler.ConflictMode.UPDATE,
        )

        self.results.itcf   = F_xt
        self.results.x_grid = wvg
        self.results.t_grid = t_grid
        self.results.run_id = run_id
        return run_id


# ============================ DB / shape helpers ============================

def _infer_theory_tag(scheme_module, rs: float, theta: float, **scheme_kwargs) -> str:
    tmp = scheme_module.Input(rs, theta, **scheme_kwargs)
    return str(getattr(tmp, "theory", scheme_module.__name__.split(".")[-1].upper()))

def _ensure_run_for(dbh: database.DataBaseHandler, scheme_module, theory: str,
                    rs: float, theta: float, desired_dim: str, **scheme_kwargs) -> int:
    rid = _find_run_id(dbh, theory, rs, theta, desired_dim)
    if rid is None:
        # Compute a fresh run with correct native enum
        kw_no_dim = {k: v for k, v in scheme_kwargs.items() if k != "dimension"}
        inputs = scheme_module.Input(rs, theta, **kw_no_dim)
        inputs.dimension = getattr(qnative.Dimension, desired_dim)  # enum, not string
        solver = scheme_module.Solver()
        solver.compute(inputs)
        if getattr(solver, "run_id", None) is None:
            raise RuntimeError("Scheme solver did not expose run_id after compute().")
        rid = int(solver.run_id)
    else:
        # sanity: confirm stored dimension matches
        inp = dbh.get_inputs(int(rid), names=["dimension"])
        run_dim = str(inp.get("dimension", "D3")).upper()
        if run_dim != desired_dim:
            raise RuntimeError(f"Run {rid} has dimension {run_dim}, expected {desired_dim}.")
    return rid

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
        run_dim = str(inp.get("dimension", "D3")).upper()
        if run_dim == desired_dim:
            candidates.append(r)
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x.get("date", ""), x.get("time", "")))
    return int(candidates[-1]["id"])

def _read_wvg_idr_lfc(dbh: database.DataBaseHandler, run_id: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    out = dbh.get_results(run_id, names=["wvg", "idr", "lfc"])
    wvg = out.get("wvg", None)
    idr = out.get("idr", None)
    lfc = out.get("lfc", None)
    if wvg is None or idr is None or lfc is None:
        present = list(dbh.get_results(run_id, names=None).keys())
        raise ValueError(f"Run {run_id} lacks 'wvg'/'idr'/'lfc' required for ITCF. Present result keys: {present}")
    return np.asarray(wvg, float), np.asarray(idr, float), np.asarray(lfc, float)

def _normalize_idr(idr: np.ndarray, nx: int) -> Tuple[np.ndarray, int]:
    """Ensure idr shape is (nx, nl). If it's (nl, nx), transpose."""
    if idr.ndim != 2:
        raise ValueError(f"Expected idr 2D, got shape {idr.shape}")
    if idr.shape[0] == nx:
        return idr, idr.shape[1]
    if idr.shape[1] == nx:
        return idr.T, idr.shape[0]
    raise ValueError(f"Cannot align idr shape {idr.shape} with wvg length {nx}")

def _normalize_lfc(lfc_arr: np.ndarray, nx: int, nl: int) -> np.ndarray:
    """
    Ensure lfc(x,l) has shape (nx, nl) by broadcasting or transposing as needed.
    Accepts (nx,), (nl,), (nx,1), (1,nx), (nl,1), (1,nl), (nx,nl), (nl,nx).
    """
    a = np.asarray(lfc_arr)
    if a.ndim == 1:
        if a.shape[0] == nx:   # lfc(x)
            return np.repeat(a[:, None], nl, axis=1)
        if a.shape[0] == nl:   # lfc(l)
            return np.repeat(a[None, :], nx, axis=0)
        raise ValueError(f"1D lfc length {a.shape[0]} matches neither nx={nx} nor nl={nl}")
    if a.ndim == 2:
        r, c = a.shape
        if (r, c) == (nx, nl): return a
        if (r, c) == (nl, nx): return a.T
        if (r, c) == (nx, 1):  return np.repeat(a, nl, axis=1)
        if (r, c) == (1, nx):  return np.repeat(a.reshape(nx)[None, :], nl, axis=0).T
        if (r, c) == (nl, 1):  return np.repeat(a, nx, axis=1).T
        if (r, c) == (1, nl):  return np.repeat(a, nx, axis=0)
        raise ValueError(f"Unsupported 2D lfc shape {a.shape} for nx={nx}, nl={nl}")
    raise ValueError(f"lfc must be 1D or 2D, got shape {a.shape}")


# =============================== 2D pieces ================================

def _chemical_potential_2d(Theta: float) -> float:
    # μ̄ = ln(e^{1/Θ} - 1)
    return np.log(np.exp(1.0 / float(Theta)) - 1.0)

def _fermi_dirac_integral_mhalf(x):
    """
    F_{-1/2}(x) = ∫_0^∞ (1/√π) t^{-1/2} / (exp(t - x) + 1) dt
    Works for scalar or array x. Uses asymptotics outside [-20, 20].
    """
    def _integrand(t, xx):
        if t <= 0.0:
            return 0.0
        return (1.0 / np.sqrt(np.pi)) * (t**(-0.5)) / (np.exp(t - xx) + 1.0)

    arr = np.asarray(x, dtype=float)
    scalar = arr.ndim == 0
    if scalar:
        arr = arr.reshape(1)

    out = np.empty_like(arr)

    # large negative: e^x * sqrt(pi)
    mask_neg = arr < -20.0
    if np.any(mask_neg):
        out[mask_neg] = np.exp(arr[mask_neg]) * np.sqrt(np.pi)

    # large positive: Sommerfeld expansion
    mask_pos = arr > 20.0
    if np.any(mask_pos):
        xp = arr[mask_pos]
        out[mask_pos] = 2.0 * np.sqrt(xp) * (
            1.0 - (np.pi**2) / (48.0 * xp) + 7.0 * (np.pi**4) / (3840.0 * xp**2)
        )

    # moderate region: quad per element (write back by absolute index!)
    mask_mid = ~(mask_neg | mask_pos)
    if np.any(mask_mid):
        idxs = np.flatnonzero(mask_mid)
        for k in idxs:
            val = float(arr[k])
            val_int, _ = integrate.quad(_integrand, 0.0, np.inf,
                                        args=(val,), limit=200, epsabs=1e-9, epsrel=1e-7)
            out[k] = val_int

    return float(out[0]) if scalar else out

def _fd_mhalf_prime(eta: np.ndarray, h: float = 1e-3) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    return (_fermi_dirac_integral_mhalf(eta + h) - _fermi_dirac_integral_mhalf(eta - h)) / (2.0 * h)

def _ratio_cosh_over_sinh_exact(z, t):
    """cosh(z*(t-1/2)) / sinh(z/2) with broadcasting."""
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    return np.cosh(z * (t - 0.5)) / np.sinh(0.5 * z)

def _leggauss_interval(n: int, ymax: float):
    """
    Nodes and weights for ∫_0^{ymax} f(y) dy ≈ dot(W, f(X)).
    """
    x, w = np.polynomial.legendre.leggauss(n)  # on [-1,1]
    X = 0.5*(x + 1.0) * ymax
    W = 0.5*ymax * w
    return X, W

def _F_HF_2d_matrix(x_grid: np.ndarray, t_grid: np.ndarray, Theta: float,
                    z_thresh: float = 5e-3,
                    n_small: int = 32,
                    n_large: int = 64,
                    y_max: float = 30.0) -> np.ndarray:
    """
    Hybrid integrator for 2D F_HF:

    ∫_0^{y_max} dy [ (F_{-1/2}(η1)-F_{-1/2}(η2)) * cosh(z(t-1/2))/sinh(z/2) ]
    with z = x*y/Θ, η1 = μ̄ - (y-x)^2/(4Θ), η2 = μ̄ - (y+x)^2/(4Θ).

    On [0, y0] (z <= z_thresh) use small-z limit:  2 * dF_{-1/2}/dη (η̄), η̄ = μ̄ - (y^2 + x^2)/(4Θ).
    On [y0, y_max] use Gauss–Legendre on the full integrand.
    Prefactor 0.5*sqrt(Θ/π) is applied to both contributions.
    """
    mu_bar = _chemical_potential_2d(Theta)
    X = np.asarray(x_grid, dtype=float)
    T = np.asarray(t_grid, dtype=float)
    nx, nt = X.size, T.size
    out = np.zeros((nx, nt), dtype=float)

    for i, x in enumerate(X):
        if x == 0.0:
            continue

        # boundary between "small-z" and "normal" regions
        y0 = min(y_max, z_thresh * Theta / max(x, 1e-300))

        # ---- small-z: t-independent approximation
        if y0 > 0.0:
            ys, ws = _leggauss_interval(n_small, y0)                   # (n_small,)
            eta_mid = mu_bar - (ys*ys + x*x) / (4.0 * Theta)           # (n_small,)
            deriv = _fd_mhalf_prime(eta_mid)                           # (n_small,)
            small_val = 2.0 * (ws @ deriv)                             # scalar
            out[i, :] += 0.5 * np.sqrt(Theta/np.pi) * small_val        # broadcast over t

        # ---- large-z: full integrand on [y0, y_max]
        if y0 < y_max:
            yl, wl = _leggauss_interval(n_large, y_max)                # (n_large,)
            mask = yl >= y0                                            # (n_large,)
            yl = yl[mask]; wl = wl[mask]
            if yl.size:
                z = (x * yl) / Theta                                   # (n_large,)
                eta1 = mu_bar - ((yl - x)**2) / (4.0 * Theta)
                eta2 = mu_bar - ((yl + x)**2) / (4.0 * Theta)
                F1 = _fermi_dirac_integral_mhalf(eta1)                 # (n_large,)
                F2 = _fermi_dirac_integral_mhalf(eta2)                 # (n_large,)
                diff  = (F1 - F2)[:, None]                              # (n_large, 1)
                ratio = _ratio_cosh_over_sinh_exact(z[:, None], T[None, :])  # (n_large, nt)
                integrand = diff * ratio
                out[i, :] += 0.5 * np.sqrt(Theta/np.pi) * (wl @ integrand)

    return out


# =============================== 3D pieces ================================

def _mu_3d(theta: float) -> float:
    theta = float(theta)
    def root_fn(mu):
        f = lambda x: np.sqrt(x) / (np.exp(x - mu) + 1.0)
        integ = integrate.quad(f, 0.0, 30.0, limit=200)[0]
        return integ - (2.0 / 3.0) * theta ** (-1.5)
    return optimize.brentq(root_fn, -100.0, 100.0, xtol=2e-5, rtol=8.881e-5, maxiter=500)

def _integrand_3d(y, x, t, theta, mu):
    if x == 0.0:
        return 0.0
    if y == 0.0:
        return 2.0 / (1.0 + np.exp((x * x) / (4.0 * theta) - mu))
    xy = x * y
    sinh_term = np.sinh(xy / (2.0 * theta))
    cosh_term = np.cosh(xy / theta * (t - 0.5))
    num = 1.0 + np.exp(mu - (x - y) ** 2 / (4.0 * theta))
    den = 1.0 + np.exp(mu - (x + y) ** 2 / (4.0 * theta))
    ln_term = np.log(num / den)
    return (cosh_term / sinh_term) * ln_term

def _F_HF_3d_matrix(x_grid: np.ndarray, t_grid: np.ndarray, theta: float, mu: float,
                    n_quad: int = 96) -> np.ndarray:
    """
    Vectorized Gauss–Legendre:
        F_HF(x,t) = (3 θ / 8) * ∫_0^{Ymax} [ln(...) * cosh/sinh] dy
    """
    X = np.asarray(x_grid, dtype=float)
    T = np.asarray(t_grid, dtype=float)
    nx, nt = X.size, T.size
    out = np.zeros((nx, nt), dtype=float)

    def y_max_for_x(x):
        base = 25.0 if not math.isclose(theta, 0.1, rel_tol=0, abs_tol=1e-12) else 13.0
        return float(np.clip(base + 8.0 / max(x, 1e-6), base, 32.0))

    pref = (3.0 * theta / 8.0)

    for i, x in enumerate(X):
        if x == 0.0:
            continue

        ymax = y_max_for_x(x)
        Y, W = _leggauss_interval(n_quad, ymax)               # (n_quad,)

        xy = x * Y
        ratio = _ratio_cosh_over_sinh_exact((xy/theta)[:, None], T[None, :])  # (n_quad, nt)

        num = 1.0 + np.exp(mu - (x - Y)**2 / (4.0 * theta))
        den = 1.0 + np.exp(mu - (x + Y)**2 / (4.0 * theta))
        ln_term = np.log(num / den)                            # (n_quad,)

        integrand = ratio * ln_term[:, None]                   # (n_quad, nt)
        out[i, :] = pref * (W @ integrand)

    return out
