from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import integrate, optimize

from qupled.schemes import hf
from qupled.util import serialize
from qupled.util.dimension import Dimension


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Inputs for ITCF postprocessing.
    Only ``coupling``, ``degeneracy``, ``dimension`` and ``dt`` are used.
    """

    dt: float = 0.01


@serialize.serializable_dataclass
class Result:
    itcf: np.ndarray = None
    x_grid: np.ndarray = None
    t_grid: np.ndarray = None


class Solver:
    """
    ITCF postprocessor acting on already-computed scheme results.
    """

    def __init__(self):
        self.results: Result = Result()

    def compute(self, inputs: Input, scheme_results: any):
        wvg = np.asarray(scheme_results.wvg, dtype=float)
        idr = np.asarray(scheme_results.idr, dtype=float)
        lfc = np.asarray(scheme_results.lfc, dtype=float)
        nx = int(wvg.shape[0])

        if inputs.dt <= 0.0:
            raise ValueError("dt must be > 0")
        t_grid = np.arange(0.0, 0.5 + inputs.dt, inputs.dt, dtype=float)
        t_grid[-1] = min(t_grid[-1], 0.5)

        idr, nl = _normalize_idr(idr, nx)
        lfc_eff = _normalize_lfc(lfc, nx, nl)

        theta = float(inputs.degeneracy)
        rs = float(inputs.coupling)
        dim = _normalize_dimension(inputs.dimension)

        if dim == "D2":
            fhf_mat = _F_HF_2d_matrix(wvg, t_grid, theta)
            norm = theta
            ip_x = np.zeros_like(wvg, dtype=float)
            mask = wvg > 0.0
            ip_x[mask] = np.sqrt(2.0) * rs / wvg[mask]
        else:
            lam3d = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
            mu3d = _mu_3d(theta)
            fhf_mat = _F_HF_3d_matrix(wvg, t_grid, theta, mu3d)
            norm = 1.5 * theta
            ip_x = np.zeros_like(wvg, dtype=float)
            mask = wvg > 0.0
            ip_x[mask] = (4.0 * lam3d * rs) / (np.pi * wvg[mask] ** 2)

        l_idx = np.arange(nl, dtype=float)
        w_l = np.ones(nl, dtype=float)
        if nl >= 1:
            w_l[1:] = 2.0
        cos_lt = np.cos(2.0 * np.pi * np.outer(l_idx, t_grid))

        one_m_lfc = 1.0 - lfc_eff
        denom = 1.0 + ip_x[:, None] * one_m_lfc * idr
        denom = np.where(np.abs(denom) < 1.0e-14, 1.0e-14, denom)
        a_xl = (idr * idr) * one_m_lfc / denom
        s_xt = (a_xl * w_l[None, :]) @ cos_lt
        f_xt = fhf_mat - (norm * ip_x[:, None]) * s_xt
        f_xt[wvg == 0.0, :] = 0.0

        self.results.itcf = f_xt
        self.results.x_grid = wvg
        self.results.t_grid = t_grid


def _normalize_dimension(dimension: Dimension | str) -> str:
    if isinstance(dimension, Dimension):
        return dimension.value
    value = str(dimension).upper()
    if value.endswith("._2D"):
        return "D2"
    if value.endswith("._3D"):
        return "D3"
    if value in {"D2", "D3"}:
        return value
    raise ValueError(f"Unsupported dimension: {dimension}")


def _normalize_idr(idr: np.ndarray, nx: int):
    if idr.ndim != 2:
        raise ValueError(f"Expected idr 2D, got shape {idr.shape}")
    if idr.shape[0] == nx:
        return idr, idr.shape[1]
    if idr.shape[1] == nx:
        return idr.T, idr.shape[0]
    raise ValueError(f"Cannot align idr shape {idr.shape} with wvg length {nx}")


def _normalize_lfc(lfc_arr: np.ndarray, nx: int, nl: int):
    a = np.asarray(lfc_arr)
    if a.ndim == 1:
        if a.shape[0] == nx:
            return np.repeat(a[:, None], nl, axis=1)
        if a.shape[0] == nl:
            return np.repeat(a[None, :], nx, axis=0)
        raise ValueError(f"1D lfc length {a.shape[0]} matches neither nx={nx} nor nl={nl}")
    if a.ndim == 2:
        r, c = a.shape
        if (r, c) == (nx, nl):
            return a
        if (r, c) == (nl, nx):
            return a.T
        if (r, c) == (nx, 1):
            return np.repeat(a, nl, axis=1)
        if (r, c) == (1, nx):
            return np.repeat(a.reshape(nx)[None, :], nl, axis=0).T
        if (r, c) == (nl, 1):
            return np.repeat(a, nx, axis=1).T
        if (r, c) == (1, nl):
            return np.repeat(a, nx, axis=0)
    raise ValueError(f"Unsupported lfc shape {a.shape} for nx={nx}, nl={nl}")


def _chemical_potential_2d(theta: float) -> float:
    return np.log(np.exp(1.0 / float(theta)) - 1.0)


def _fermi_dirac_integral_mhalf(x):
    def _integrand(t, xx):
        if t <= 0.0:
            return 0.0
        return (1.0 / np.sqrt(np.pi)) * (t ** (-0.5)) / (np.exp(t - xx) + 1.0)

    arr = np.asarray(x, dtype=float)
    scalar = arr.ndim == 0
    if scalar:
        arr = arr.reshape(1)
    out = np.empty_like(arr)
    mask_neg = arr < -20.0
    if np.any(mask_neg):
        out[mask_neg] = np.exp(arr[mask_neg]) * np.sqrt(np.pi)
    mask_pos = arr > 20.0
    if np.any(mask_pos):
        xp = arr[mask_pos]
        out[mask_pos] = 2.0 * np.sqrt(xp) * (
            1.0 - (np.pi**2) / (48.0 * xp) + 7.0 * (np.pi**4) / (3840.0 * xp**2)
        )
    mask_mid = ~(mask_neg | mask_pos)
    if np.any(mask_mid):
        idxs = np.flatnonzero(mask_mid)
        for k in idxs:
            val = float(arr[k])
            val_int, _ = integrate.quad(
                _integrand, 0.0, np.inf, args=(val,), limit=200, epsabs=1.0e-9, epsrel=1.0e-7
            )
            out[k] = val_int
    return float(out[0]) if scalar else out


def _fd_mhalf_prime(eta: np.ndarray, h: float = 1.0e-3) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    return (_fermi_dirac_integral_mhalf(eta + h) - _fermi_dirac_integral_mhalf(eta - h)) / (2.0 * h)


def _ratio_cosh_over_sinh_exact(z, t):
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    return np.cosh(z * (t - 0.5)) / np.sinh(0.5 * z)


def _leggauss_interval(n: int, ymax: float):
    x, w = np.polynomial.legendre.leggauss(n)
    X = 0.5 * (x + 1.0) * ymax
    W = 0.5 * ymax * w
    return X, W


def _F_HF_2d_matrix(
    x_grid: np.ndarray,
    t_grid: np.ndarray,
    theta: float,
    z_thresh: float = 5.0e-3,
    n_small: int = 32,
    n_large: int = 64,
    y_max: float = 30.0,
):
    mu_bar = _chemical_potential_2d(theta)
    X = np.asarray(x_grid, dtype=float)
    T = np.asarray(t_grid, dtype=float)
    nx, nt = X.size, T.size
    out = np.zeros((nx, nt), dtype=float)

    for i, x in enumerate(X):
        if x == 0.0:
            continue
        y0 = min(y_max, z_thresh * theta / max(x, 1.0e-300))
        if y0 > 0.0:
            ys, ws = _leggauss_interval(n_small, y0)
            eta_mid = mu_bar - (ys * ys + x * x) / (4.0 * theta)
            deriv = _fd_mhalf_prime(eta_mid)
            small_val = 2.0 * (ws @ deriv)
            out[i, :] += 0.5 * np.sqrt(theta / np.pi) * small_val
        if y0 < y_max:
            yl, wl = _leggauss_interval(n_large, y_max)
            mask = yl >= y0
            yl = yl[mask]
            wl = wl[mask]
            if yl.size:
                z = (x * yl) / theta
                eta1 = mu_bar - ((yl - x) ** 2) / (4.0 * theta)
                eta2 = mu_bar - ((yl + x) ** 2) / (4.0 * theta)
                F1 = _fermi_dirac_integral_mhalf(eta1)
                F2 = _fermi_dirac_integral_mhalf(eta2)
                diff = (F1 - F2)[:, None]
                ratio = _ratio_cosh_over_sinh_exact(z[:, None], T[None, :])
                out[i, :] += 0.5 * np.sqrt(theta / np.pi) * (wl @ (diff * ratio))
    return out


def _mu_3d(theta: float) -> float:
    theta = float(theta)

    def root_fn(mu):
        f = lambda x: np.sqrt(x) / (np.exp(x - mu) + 1.0)
        integ = integrate.quad(f, 0.0, 30.0, limit=200)[0]
        return integ - (2.0 / 3.0) * theta ** (-1.5)

    return optimize.brentq(root_fn, -100.0, 100.0, xtol=2.0e-5, rtol=8.881e-5, maxiter=500)


def _F_HF_3d_matrix(x_grid: np.ndarray, t_grid: np.ndarray, theta: float, mu: float, n_quad: int = 96):
    X = np.asarray(x_grid, dtype=float)
    T = np.asarray(t_grid, dtype=float)
    nx, nt = X.size, T.size
    out = np.zeros((nx, nt), dtype=float)

    def y_max_for_x(x):
        base = 25.0 if not math.isclose(theta, 0.1, rel_tol=0.0, abs_tol=1.0e-12) else 13.0
        return float(np.clip(base + 8.0 / max(x, 1.0e-6), base, 32.0))

    pref = 3.0 * theta / 8.0
    for i, x in enumerate(X):
        if x == 0.0:
            continue
        ymax = y_max_for_x(x)
        Y, W = _leggauss_interval(n_quad, ymax)
        xy = x * Y
        ratio = _ratio_cosh_over_sinh_exact((xy / theta)[:, None], T[None, :])
        num = 1.0 + np.exp(mu - (x - Y) ** 2 / (4.0 * theta))
        den = 1.0 + np.exp(mu - (x + Y) ** 2 / (4.0 * theta))
        ln_term = np.log(num / den)
        out[i, :] = pref * (W @ (ratio * ln_term[:, None]))
    return out
