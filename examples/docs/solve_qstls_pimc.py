import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

from qupled.schemes import esa, qstls, qstlsf0, qstlspimc


LAMBDA = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)


def rs_label(rs: float) -> str:
    if np.isclose(rs, round(rs), atol=1.0e-12):
        return str(int(round(rs)))
    return str(rs).replace(".", "p")


def rs_display(rs: float) -> str:
    return f"{rs:g}"


def resolve_state_dir(pimc_root: Path, rs: float) -> Path | None:
    candidates = [
        pimc_root / f"N34_rs{rs_label(rs)}_theta1",
        pimc_root / f"N34_rs{rs:g}_theta1",
        pimc_root / f"N34_rs{int(round(rs))}_theta1",
        pimc_root / f"N14_rs{rs_label(rs)}_theta1",
        pimc_root / f"N14_rs{rs:g}_theta1",
        pimc_root / f"N14_rs{int(round(rs))}_theta1",
    ]
    for path in candidates:
        if path.exists():
            return path
    pattern = f"*rs{rs:g}_theta1"
    matches = sorted(pimc_root.glob(pattern))
    if matches:
        return matches[0]
    return None


def run_state_point(rs: float, pimc_eta: float):
    scheme = qstlspimc.Solver()
    inputs = qstlspimc.Input(
        coupling=rs,
        degeneracy=1.0,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
        pimc_eta=pimc_eta,
    )
    scheme.compute(inputs)
    print(
        f"QSTLS-PIMC: rs={rs_display(rs)}, eta={pimc_eta:.3g}, uint={scheme.results.uint:.8f}"
    )
    f0_grid = getattr(scheme.results, "f0_grid", None)
    f0_values = getattr(scheme.results, "f0_values", None)
    return (
        scheme.results.wvg,
        scheme.results.ssf,
        scheme.results.lfc[:, 0],
        f0_grid,
        f0_values,
    )


def run_f0_state_point(rs: float, pimc_eta: float):
    scheme = qstlsf0.Solver()
    inputs = qstlsf0.Input(
        coupling=rs,
        degeneracy=1.0,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
        pimc_eta=pimc_eta,
    )
    scheme.compute(inputs)
    print(
        f"QSTLS-F0: rs={rs_display(rs)}, eta={pimc_eta:.3g}, uint={scheme.results.uint:.8f}"
    )
    f0_grid = getattr(scheme.results, "f0_grid", None)
    f0_values = getattr(scheme.results, "f0_values", None)
    return (
        scheme.results.wvg,
        scheme.results.ssf,
        scheme.results.lfc[:, 0],
        f0_grid,
        f0_values,
    )


def run_esa_state_point(rs: float):
    scheme = esa.Solver()
    inputs = esa.Input(
        coupling=rs,
        degeneracy=1.0,
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"ESA: rs={rs_display(rs)}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf, scheme.results.lfc[:, 0]


def run_qstls_state_point(rs: float):
    scheme = qstls.Solver()
    inputs = qstls.Input(
        coupling=rs,
        degeneracy=1.0,
        matsubara=300,
        mixing=0.5,
        threads=16,
        integral_strategy="segregated",
        cutoff=10.0,
        resolution=0.2,
    )
    scheme.compute(inputs)
    print(f"QSTLS: rs={rs_display(rs)}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf, scheme.results.lfc[:, 0]


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
    state_dir = resolve_state_dir(pimc_root, rs)
    if state_dir is None:
        print(f"Missing PIMC state directory for rs={rs_display(rs)} under {pimc_root}")
        return None
    key_path = state_dir / "static_structure_factor_key.dat"
    ssf_values = []
    if key_path.exists():
        key_data = np.loadtxt(key_path)
        q_index = key_data[:, 0].astype(int)
        q_values = key_data[:, 1]
        for idx in q_index:
            tau_path = state_dir / f"Tau_k_index{idx}fermion_density_response_electron_34.res"
            if not tau_path.exists():
                print(f"Missing PIMC tau file: {tau_path}")
                return None
            tau_data = np.loadtxt(tau_path)
            tau_zero = np.atleast_2d(tau_data)[0, 1]
            ssf_values.append(tau_zero)
    else:
        intfive_files = sorted(state_dir.glob("IntFive_fermion_density_response_electron_*.res"))
        tau_files = sorted(
            state_dir.glob("Tau_k_index*fermion_density_response_electron_*.res"),
            key=lambda path: int(path.name.split("Tau_k_index", 1)[1].split("fermion", 1)[0]),
        )
        if not intfive_files or not tau_files:
            print(f"Missing PIMC key/intfive data in {state_dir}")
            return None
        intfive_data = np.loadtxt(intfive_files[0])
        q_values = np.atleast_2d(intfive_data)[:, 0]
        npts = min(len(q_values), len(tau_files))
        q_values = q_values[:npts]
        for tau_path in tau_files[:npts]:
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
    state_dir = resolve_state_dir(pimc_root, rs)
    if state_dir is None:
        print(f"Missing PIMC state directory for rs={rs_display(rs)} under {pimc_root}")
        return None
    vector_path = state_dir / "Vector_iFreq0.res"
    if not vector_path.exists():
        print(f"Missing PIMC vector file: {vector_path}")
        return None

    data = np.loadtxt(vector_path)
    q_values = data[:, 0]
    if data.shape[1] > 9:
        lfc_values = data[:, 9]
    elif data.shape[1] >= 4:
        lfc_values = data[:, 3]
    else:
        print(f"Unsupported PIMC vector layout in {vector_path}: shape={data.shape}")
        return None
    x_values = LAMBDA * rs * q_values
    return x_values, lfc_values


def load_pimc_momentum(rs: float):
    path = Path.cwd() / f"Momentum_UEG_rs{rs_label(rs)}_theta1.txt"
    if not path.exists():
        alt = Path("/home/fotios/qupled") / path.name
        path = alt if alt.exists() else path
    if not path.exists():
        print(f"Missing PIMC momentum file: {path}")
        return None

    rows = []
    rows_i = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            x = float(parts[0])
            y = max(float(parts[1]), 0.0)
            rows.append((x, y))
            if len(parts) >= 5 and parts[4].lower() == "i":
                rows_i.append((x, y))
    data = np.array(rows_i if len(rows_i) >= 4 else rows, dtype=float)
    x_values = LAMBDA * rs * data[:, 0]
    return x_values, data[:, 1]


def fermi_dirac_half_integral(mu: float) -> float:
    t_max = max(60.0, mu + 60.0)
    t = np.linspace(0.0, t_max, 12000)
    integrand = np.sqrt(t) / (np.exp(t - mu) + 1.0)
    return np.trapz(integrand, t)


def chemical_potential_3d(theta: float = 1.0) -> float:
    target = 2.0 / (3.0 * theta**1.5)
    a, b = -20.0, 20.0
    fa = fermi_dirac_half_integral(a) - target
    fb = fermi_dirac_half_integral(b) - target
    while fa * fb > 0.0:
        a -= 10.0
        b += 10.0
        fa = fermi_dirac_half_integral(a) - target
        fb = fermi_dirac_half_integral(b) - target
    for _ in range(80):
        c = 0.5 * (a + b)
        fc = fermi_dirac_half_integral(c) - target
        if fa * fc <= 0.0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


def fermi_dirac_distribution(y: np.ndarray, theta: float = 1.0) -> np.ndarray:
    mu = chemical_potential_3d(theta)
    return 1.0 / (np.exp(y * y / theta - mu) + 1.0)


def main():
    rs_values = (3.23, 10.0, 20.0)
    eta_values = (0.5, 1.0, 2.0)
    qstls_pimc_data = {
        rs: {eta: run_state_point(rs, eta) for eta in eta_values} for rs in rs_values
    }
    qstls_f0_data = {
        rs: {eta: run_f0_state_point(rs, eta) for eta in eta_values} for rs in rs_values
    }
    qstls_data = {rs: run_qstls_state_point(rs) for rs in rs_values}
    esa_data = {rs: run_esa_state_point(rs) for rs in rs_values}
    pimc_data = {rs: load_pimc_ssf(rs) for rs in rs_values}
    pimc_lfc_data = {rs: load_pimc_lfc(rs) for rs in rs_values}
    pimc_momentum_data = {rs: load_pimc_momentum(rs) for rs in rs_values}

    for rs in rs_values:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.2))

        for eta, (wvg, ssf, lfc0, f0_grid, f0_values) in qstls_pimc_data[rs].items():
            ax1.plot(wvg, ssf, label=fr"QSTLS-PIMC $\eta={eta:g}$")
            ax2.plot(wvg, lfc0, label=fr"QSTLS-PIMC $\eta={eta:g}$")
            if (
                f0_grid is not None
                and f0_values is not None
                and len(f0_grid)
                and len(f0_values)
            ):
                ax3.plot(f0_grid, f0_values, label=fr"QSTLS-PIMC $\eta={eta:g}$")
            else:
                print(
                    f"Missing native f0 output for rs={rs_display(rs)}, eta={eta:g}. "
                    "Rebuild the editable install to expose f0_grid/f0_values."
                )
        for eta, (wvg, ssf, lfc0, f0_grid, f0_values) in qstls_f0_data[rs].items():
            ax1.plot(wvg, ssf, "--", label=fr"QSTLS-F0 $\eta={eta:g}$")
            ax2.plot(wvg, lfc0, "--", label=fr"QSTLS-F0 $\eta={eta:g}$")
            if (
                f0_grid is not None
                and f0_values is not None
                and len(f0_grid)
                and len(f0_values)
            ):
                ax3.plot(f0_grid, f0_values, "--", label=fr"QSTLS-F0 $\eta={eta:g}$")
            else:
                print(
                    f"Missing native f0 output for QSTLS-F0 rs={rs_display(rs)}, eta={eta:g}. "
                    "Rebuild the editable install to expose f0_grid/f0_values."
                )

        wvg_qstls, ssf_qstls, lfc_qstls = qstls_data[rs]
        ax1.plot(wvg_qstls, ssf_qstls, ":", label="QSTLS")
        ax2.plot(wvg_qstls, lfc_qstls, ":", label="QSTLS")

        wvg_esa, ssf_esa, lfc_esa = esa_data[rs]
        ax1.plot(wvg_esa, ssf_esa, "-.", label="ESA")
        ax2.plot(wvg_esa, lfc_esa, "-.", label="ESA")

        pimc_ssf = pimc_data[rs]
        if pimc_ssf is not None:
            x, y = pimc_ssf
            ax1.scatter(x, y, s=22, alpha=0.75, label="PIMC")

        pimc_lfc = pimc_lfc_data[rs]
        if pimc_lfc is not None:
            x, y = pimc_lfc
            ax2.scatter(x, y, s=22, alpha=0.75, label="PIMC")

        pimc_momentum = pimc_momentum_data[rs]
        if pimc_momentum is not None:
            x, y = pimc_momentum
            ax3.scatter(x, y, s=16, alpha=0.55, label="PIMC")
            y_plot = np.linspace(0.0, max(8.0, 1.1 * np.max(x)), 800)
            ax3.plot(
                y_plot,
                fermi_dirac_distribution(y_plot, 1.0),
                "--",
                lw=1.5,
                label="Fermi-Dirac",
            )

        ax1.set_xlabel("x")
        ax1.set_ylabel("S(x)")
        ax1.set_title(fr"Static structure factor, $r_s={rs_display(rs)}$")
        ax1.grid(alpha=0.3)
        ax1.legend()

        ax2.set_xlabel("x")
        ax2.set_ylabel("G(x, l=0)")
        ax2.set_title(fr"Local-field correction, $r_s={rs_display(rs)}$")
        ax2.grid(alpha=0.3)
        ax2.legend()

        ax3.set_xlabel("y")
        ax3.set_ylabel(r"$f_0(y)$")
        ax3.set_title(fr"Momentum distribution, $r_s={rs_display(rs)}$")
        ax3.grid(alpha=0.3)
        ax3.set_yscale("log")
        ax3.legend()

        fig.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
