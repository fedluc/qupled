import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

from qupled import itcf
from qupled.schemes import esa, qstls, qstlspimc


LAMBDA = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
RS = 10.0
THETA = 1.0


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


def load_pimc_itcf_slice(rs: float, tau_reduced: float):
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

    # The supplied PIMC tau files span [0, 2 beta] and are symmetric around beta.
    # Match the reduced ITCF time t in [0, 1/2] via tau_file = 2 beta t.
    beta_ha = (LAMBDA * rs) ** 2 / THETA
    tau_target = 2.0 * beta_ha * tau_reduced
    values = []

    for idx in q_index:
        tau_path = state_dir / f"Tau_k_index{idx}fermion_density_response_electron_34.res"
        if not tau_path.exists():
            print(f"Missing PIMC tau file: {tau_path}")
            return None
        tau_data = np.loadtxt(tau_path)
        tau_data = np.atleast_2d(tau_data)
        row = int(np.argmin(np.abs(tau_data[:, 0] - tau_target)))
        values.append(tau_data[row, 1])

    x_values = LAMBDA * rs * q_values
    return x_values, np.asarray(values)


def run_qstls_pimc(rs: float):
    scheme = qstlspimc.Solver()
    inputs = qstlspimc.Input(
        coupling=rs,
        degeneracy=THETA,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"QSTLS-PIMC: rs={rs:.0f}, uint={scheme.results.uint:.8f}")
    return inputs, scheme.results


def run_qstls(rs: float):
    scheme = qstls.Solver()
    inputs = qstls.Input(
        coupling=rs,
        degeneracy=THETA,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"QSTLS: rs={rs:.0f}, uint={scheme.results.uint:.8f}")
    return inputs, scheme.results


def run_esa(rs: float):
    scheme = esa.Solver()
    inputs = esa.Input(
        coupling=rs,
        degeneracy=THETA,
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"ESA: rs={rs:.0f}, uint={scheme.results.uint:.8f}")
    return inputs, scheme.results


def compute_itcf(rs: float, scheme_inputs, scheme_results):
    solver = itcf.Solver()
    inputs = itcf.Input(
        coupling=rs,
        degeneracy=THETA,
        dimension=scheme_inputs.dimension,
        dt=0.01,
    )
    solver.compute(inputs, scheme_results)
    return solver.results


def get_itcf_slice(itcf_results, tau_reduced: float):
    idx = int(np.argmin(np.abs(itcf_results.t_grid - tau_reduced)))
    return itcf_results.x_grid, itcf_results.itcf[:, idx]


def main():
    qstls_pimc_inputs, qstls_pimc_results = run_qstls_pimc(RS)
    qstls_inputs, qstls_results = run_qstls(RS)
    esa_inputs, esa_results = run_esa(RS)

    qstls_pimc_itcf = compute_itcf(RS, qstls_pimc_inputs, qstls_pimc_results)
    qstls_itcf = compute_itcf(RS, qstls_inputs, qstls_results)
    esa_itcf = compute_itcf(RS, esa_inputs, esa_results)

    pimc_ssf = load_pimc_itcf_slice(RS, 0.0)
    pimc_tau_quarter = load_pimc_itcf_slice(RS, 0.25)
    pimc_tau_half = load_pimc_itcf_slice(RS, 0.5)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.2))

    ax1.plot(qstls_pimc_results.wvg, qstls_pimc_results.ssf, label="QSTLS-PIMC")
    ax1.plot(qstls_results.wvg, qstls_results.ssf, "--", label="QSTLS")
    ax1.plot(esa_results.wvg, esa_results.ssf, ":", label="ESA")
    if pimc_ssf is not None:
        x, y = pimc_ssf
        ax1.scatter(x, y, s=24, alpha=0.75, label="PIMC")

    x, y = get_itcf_slice(qstls_pimc_itcf, 0.25)
    ax2.plot(x, y, label="QSTLS-PIMC")
    x, y = get_itcf_slice(qstls_itcf, 0.25)
    ax2.plot(x, y, "--", label="QSTLS")
    x, y = get_itcf_slice(esa_itcf, 0.25)
    ax2.plot(x, y, ":", label="ESA")
    if pimc_tau_quarter is not None:
        x, y = pimc_tau_quarter
        ax2.scatter(x, y, s=24, alpha=0.75, label="PIMC")

    x, y = get_itcf_slice(qstls_pimc_itcf, 0.5)
    ax3.plot(x, y, label="QSTLS-PIMC")
    x, y = get_itcf_slice(qstls_itcf, 0.5)
    ax3.plot(x, y, "--", label="QSTLS")
    x, y = get_itcf_slice(esa_itcf, 0.5)
    ax3.plot(x, y, ":", label="ESA")
    if pimc_tau_half is not None:
        x, y = pimc_tau_half
        ax3.scatter(x, y, s=24, alpha=0.75, label="PIMC")

    ax1.set_xlabel("x")
    ax1.set_ylabel("S(x)")
    ax1.set_title(r"Static Structure Factor, $r_s=10$")
    ax1.grid(alpha=0.3)
    ax1.legend()

    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$F(x,\tau)$")
    ax2.set_title(r"ITCF, $\tau=1/4$")
    ax2.grid(alpha=0.3)
    ax2.legend()

    ax3.set_xlabel("x")
    ax3.set_ylabel(r"$F(x,\tau)$")
    ax3.set_title(r"ITCF, $\tau=1/2$")
    ax3.grid(alpha=0.3)
    ax3.legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
