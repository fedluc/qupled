import matplotlib.pyplot as plt

from qupled.database.base_tables import RunStatus, TableKeys as BaseTableKeys
from qupled.database.scheme_tables import TableKeys
from qupled.postprocess.output import DataBase, OutputType
from qupled.schemes import esa, qstls, ve


RS = 10.0
THETA = 1.0
MATSUBARA = 64
MIXING = 0.5
THREADS = 16
RESOLUTION = 0.1
CUTOFF = 10.0
COUPLING_RESOLUTION = 0.1


def _same_float(a: float, b: float, atol: float = 1.0e-14) -> bool:
    return abs(float(a) - float(b)) <= atol


def _load_matching_run(
    theory: str,
    rs: float,
    theta: float,
    required_inputs: dict,
):
    runs = DataBase.inspect_runs(OutputType.SCHEME)
    for run in runs:
        if run[BaseTableKeys.STATUS.value] != RunStatus.SUCCESS.value:
            continue
        if run[TableKeys.THEORY.value] != theory:
            continue
        if not _same_float(run[TableKeys.COUPLING.value], rs):
            continue
        if not _same_float(run[TableKeys.DEGENERACY.value], theta):
            continue
        run_id = run[BaseTableKeys.PRIMARY_KEY.value]
        run_inputs = DataBase.read_inputs(OutputType.SCHEME, run_id)
        same_inputs = True
        for key, expected in required_inputs.items():
            actual = run_inputs.get(key)
            if isinstance(expected, float):
                if actual is None or not _same_float(actual, expected):
                    same_inputs = False
                    break
            else:
                if actual != expected:
                    same_inputs = False
                    break
        if not same_inputs:
            continue
        results = DataBase.read_results(
            OutputType.SCHEME, run_id, names=["wvg", "ssf", "lfc", "uint"]
        )
        print(f"Loading {theory} from database for run_id = {run_id}")
        return results["wvg"], results["ssf"], results["lfc"][:, 0], results["uint"]
    return None


def run_qstls(rs: float, theta: float):
    cached = _load_matching_run(
        "QSTLS",
        rs,
        theta,
        {
            "matsubara": MATSUBARA,
            "mixing": MIXING,
            "cutoff": CUTOFF,
            "resolution": RESOLUTION,
            "integral_strategy": "segregated",
        },
    )
    if cached is not None:
        wvg, ssf, lfc0, uint = cached
        print(f"QSTLS: rs={rs:g}, theta={theta:g}, uint={uint:.8f}")
        return wvg, ssf, lfc0

    scheme = qstls.Solver()
    inputs = qstls.Input(
        coupling=rs,
        degeneracy=theta,
        mixing=MIXING,
        matsubara=MATSUBARA,
        threads=THREADS,
        integral_strategy="segregated",
        cutoff=CUTOFF,
        resolution=RESOLUTION,
    )
    scheme.compute(inputs)
    print(f"QSTLS: rs={rs:g}, theta={theta:g}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf, scheme.results.lfc[:, 0]


def run_ve(rs: float, theta: float):
    cached = _load_matching_run(
        "VE",
        rs,
        theta,
        {
            "matsubara": MATSUBARA,
            "mixing": MIXING,
            "cutoff": CUTOFF,
            "resolution": RESOLUTION,
            "integral_strategy": "segregated",
            "coupling_resolution": COUPLING_RESOLUTION,
        },
    )
    if cached is not None:
        wvg, ssf, lfc0, uint = cached
        print(f"VE: rs={rs:g}, theta={theta:g}, uint={uint:.8f}")
        raise RuntimeError(
            "Cached VE runs do not include diagnostics. "
            "Delete the old VE run or force a fresh solve for diagnostics."
        )

    scheme = ve.Solver()
    inputs = ve.Input(
        coupling=rs,
        degeneracy=theta,
        mixing=MIXING,
        matsubara=MATSUBARA,
        threads=THREADS,
        integral_strategy="segregated",
        cutoff=CUTOFF,
        resolution=RESOLUTION,
        coupling_resolution=COUPLING_RESOLUTION,
    )
    scheme.compute(inputs)
    print(f"VE: rs={rs:g}, theta={theta:g}, uint={scheme.results.uint:.8f}")
    return (
        scheme.results.wvg,
        scheme.results.ssf,
        scheme.results.lfc[:, 0],
        scheme.results,
    )


def run_esa(rs: float, theta: float):
    cached = _load_matching_run(
        "ESA",
        rs,
        theta,
        {
            "cutoff": CUTOFF,
            "resolution": RESOLUTION,
        },
    )
    if cached is not None:
        wvg, ssf, _, uint = cached
        print(f"ESA: rs={rs:g}, theta={theta:g}, uint={uint:.8f}")
        return wvg, ssf

    scheme = esa.Solver()
    inputs = esa.Input(
        coupling=rs,
        degeneracy=theta,
        cutoff=CUTOFF,
        resolution=RESOLUTION,
    )
    scheme.compute(inputs)
    print(f"ESA: rs={rs:g}, theta={theta:g}, uint={scheme.results.uint:.8f}")
    return scheme.results.wvg, scheme.results.ssf


def main():
    wvg_qstls, ssf_qstls, lfc_qstls = run_qstls(RS, THETA)
    wvg_ve, ssf_ve, lfc_ve, ve_results = run_ve(RS, THETA)
    wvg_esa, ssf_esa = run_esa(RS, THETA)

    print("VE diagnostics:")
    print(f"  a0      = {ve_results.a0_coeff:.8e}")
    print(f"  Kxc0    = {ve_results.kxc0:.8e}")
    print(f"  Mxc_inf = {ve_results.mxc_inf:.8e}")
    print(f"  Cxc     = {ve_results.cxc:.8e}")
    print(f"  OmegaM2 = {ve_results.omega_m2:.8e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax1.plot(wvg_qstls, ssf_qstls, "--", lw=2.0, label="QSTLS")
    ax1.plot(wvg_ve, ssf_ve, lw=2.0, label="VE")
    ax1.plot(wvg_esa, ssf_esa, ":", lw=2.0, label="ESA")
    ax1.set_xlabel("x")
    ax1.set_ylabel("S(x)")
    ax1.set_title(fr"Static structure factor, $r_s={RS:g}$, $\Theta={THETA:g}$")
    ax1.grid(alpha=0.3)
    ax1.legend()

    ax2.plot(wvg_qstls, lfc_qstls, "--", lw=2.0, label="QSTLS")
    ax2.plot(wvg_ve, lfc_ve, lw=2.0, label=r"VE, $l=0$")
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$G(x, l=0)$")
    ax2.set_title(fr"Local-field correction, $r_s={RS:g}$, $\Theta={THETA:g}$")
    ax2.grid(alpha=0.3)
    ax2.legend()

    fig.tight_layout()

    fig_diag, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))
    mats = ve_results.matsubara_grid

    axes[0].plot(mats, ve_results.a_coeff, lw=2.0)
    axes[0].set_xlabel(r"$l$")
    axes[0].set_ylabel(r"$a_l$")
    axes[0].set_title("VE coefficient")
    axes[0].grid(alpha=0.3)

    axes[1].plot(mats, ve_results.a1_coeff, lw=2.0)
    axes[1].set_xlabel(r"$l$")
    axes[1].set_ylabel(r"$A_1(l)$")
    axes[1].set_title("Small-x kernel coefficient")
    axes[1].grid(alpha=0.3)

    axes[2].plot(mats, ve_results.mxc_l, lw=2.0)
    axes[2].axhline(ve_results.kxc0, ls="--", color="tab:red", label=r"$K_{xc}^0$")
    axes[2].axhline(
        ve_results.mxc_inf, ls=":", color="tab:green", label=r"$M_{xc}^{\infty}$"
    )
    axes[2].set_xlabel(r"$l$")
    axes[2].set_ylabel(r"$M_{xc}(l)$")
    axes[2].set_title("Viscoelastic modulus")
    axes[2].grid(alpha=0.3)
    axes[2].legend()

    fig_diag.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
