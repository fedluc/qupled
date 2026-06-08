import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

from qupled.schemes import esa, qvsstls, qvsstlsf0


LAMBDA = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)


def run_qvsstls_f0_state_point(rs: float, eta: float):
    scheme = qvsstlsf0.Solver()
    inputs = qvsstlsf0.Input(
        coupling=rs,
        degeneracy=1.0,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
        coupling_resolution=0.1,
        degeneracy_resolution=0.1,
        pimc_eta=eta,
    )
    scheme.compute(inputs)
    print(f"QVSSTLS-F0: rs={rs:.0f}, eta={eta:.3g}, uint={scheme.results.uint:.8f}")
    f0_grid = getattr(scheme.results, "f0_grid", None)
    f0_values = getattr(scheme.results, "f0_values", None)
    return (
        scheme.results.wvg,
        scheme.results.ssf,
        scheme.results.lfc[:, 0],
        f0_grid,
        f0_values,
    )


def run_qvsstls_state_point(rs: float):
    scheme = qvsstls.Solver()
    inputs = qvsstls.Input(
        coupling=rs,
        degeneracy=1.0,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
        coupling_resolution=0.1,
        degeneracy_resolution=0.1,
    )
    scheme.compute(inputs)
    print(f"QVSSTLS: rs={rs:.0f}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf, scheme.results.lfc[:, 0]


def run_esa_state_point(rs: float):
    scheme = esa.Solver()
    inputs = esa.Input(
        coupling=rs,
        degeneracy=1.0,
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"ESA: rs={rs:.0f}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf


def find_pimc_root():
    candidates = [
        Path.cwd() / "PIMCforqSTLS",
        Path.cwd().parent / "PIMCforqSTLS",
        Path("/home/fotios/qSTLSPIMC/PIMCforqSTLS"),
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def load_pimc_ssf(rs: float):
    pimc_root = find_pimc_root()
    if pimc_root is None:
        print("Missing PIMCforqSTLS directory.")
        return None
    state_dir = pimc_root / f"N34_rs{int(rs)}_theta1"
    key_path = state_dir / "static_structure_factor_key.dat"
    if not key_path.exists():
        print(f"Missing PIMC key file: {key_path}")
        return None

    key_data = np.loadtxt(key_path)
    q_index = key_data[:, 0].astype(int)
    q_values = key_data[:, 1]
    ssf_values = []

    for idx in q_index:
        tau_path = state_dir / f"Tau_k_index{idx}fermion_density_response_electron_34.res"
        if not tau_path.exists():
            print(f"Missing PIMC tau file: {tau_path}")
            return None
        tau_data = np.loadtxt(tau_path)
        tau_zero = np.atleast_2d(tau_data)[0, 1]
        ssf_values.append(tau_zero)

    x_values = LAMBDA * rs * q_values
    return x_values, np.asarray(ssf_values)


def load_pimc_lfc(rs: float):
    pimc_root = find_pimc_root()
    if pimc_root is None:
        print("Missing PIMCforqSTLS directory.")
        return None
    state_dir = pimc_root / f"N34_rs{int(rs)}_theta1"
    vector_path = state_dir / "Vector_iFreq0.res"
    if not vector_path.exists():
        print(f"Missing PIMC vector file: {vector_path}")
        return None

    data = np.loadtxt(vector_path)
    q_values = data[:, 0]
    lfc_values = data[:, 9]
    x_values = LAMBDA * rs * q_values
    return x_values, lfc_values


def main():
    rs_values = (10.0, 20.0)
    eta_values = (0.5, 1.0, 2.0, 4.0)

    qvsstls_f0_data = {
        rs: {eta: run_qvsstls_f0_state_point(rs, eta) for eta in eta_values}
        for rs in rs_values
    }
    qvsstls_data = {rs: run_qvsstls_state_point(rs) for rs in rs_values}
    esa_data = {rs: run_esa_state_point(rs) for rs in rs_values}
    pimc_ssf_data = {rs: load_pimc_ssf(rs) for rs in rs_values}
    pimc_lfc_data = {rs: load_pimc_lfc(rs) for rs in rs_values}

    for rs in rs_values:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.2))

        for eta, (wvg, ssf, lfc0, f0_grid, f0_values) in qvsstls_f0_data[rs].items():
            ax1.plot(wvg, ssf, label=fr"QVSSTLS-F0 $\\eta={eta:g}$")
            ax2.plot(wvg, lfc0, label=fr"QVSSTLS-F0 $\\eta={eta:g}$")
            if (
                f0_grid is not None
                and f0_values is not None
                and len(f0_grid)
                and len(f0_values)
            ):
                ax3.plot(f0_grid, f0_values, label=fr"QVSSTLS-F0 $\\eta={eta:g}$")

        wvg_qvs, ssf_qvs, lfc_qvs = qvsstls_data[rs]
        ax1.plot(wvg_qvs, ssf_qvs, "--", label="QVSSTLS")
        ax2.plot(wvg_qvs, lfc_qvs, "--", label="QVSSTLS")

        wvg_esa, ssf_esa = esa_data[rs]
        ax1.plot(wvg_esa, ssf_esa, ":", label="ESA")

        pimc_ssf = pimc_ssf_data[rs]
        if pimc_ssf is not None:
            x, y = pimc_ssf
            ax1.scatter(x, y, s=22, alpha=0.75, label="PIMC")

        pimc_lfc = pimc_lfc_data[rs]
        if pimc_lfc is not None:
            x, y = pimc_lfc
            ax2.scatter(x, y, s=22, alpha=0.75, label="PIMC")

        ax1.set_xlabel("x")
        ax1.set_ylabel("S(x)")
        ax1.set_title(fr"Static structure factor, $r_s={int(rs)}$")
        ax1.grid(alpha=0.3)
        ax1.legend()

        ax2.set_xlabel("x")
        ax2.set_ylabel("G(x, l=0)")
        ax2.set_title(fr"Local-field correction, $r_s={int(rs)}$")
        ax2.grid(alpha=0.3)
        ax2.legend()

        ax3.set_xlabel("y")
        ax3.set_ylabel(r"$f_0(y)$")
        ax3.set_title(fr"Momentum distribution from native output, $r_s={int(rs)}$")
        ax3.grid(alpha=0.3)
        ax3.set_yscale("log")
        ax3.legend()

        fig.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
