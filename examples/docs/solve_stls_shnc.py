import matplotlib.pyplot as plt

from qupled.schemes import stls, stlsiet
from qupled.util.dimension import Dimension


def run_stls(rs: float, theta: float):
    scheme = stls.Solver()
    inputs = stls.Input(
        coupling=rs,
        degeneracy=theta,
        dimension=Dimension._2D,
        mixing=0.1,
        theory="STLS",
    )
    scheme.compute(inputs)
    return scheme.results.wvg, scheme.results.ssf, scheme.results.lfc[:, 0]


def run_stls_iet(rs: float, theta: float, theory: str):
    scheme = stlsiet.Solver()
    inputs = stlsiet.Input(
        coupling=rs,
        degeneracy=theta,
        dimension=Dimension._2D,
        mixing=0.1,
        theory=theory,
    )
    scheme.compute(inputs)
    return (
        scheme.results.wvg,
        scheme.results.ssf,
        scheme.results.lfc[:, 0],
        scheme.results.bf,
    )


def main():
    rs = 10.0
    theta = 1.0

    wvg_stls, ssf_stls, lfc_stls = run_stls(rs, theta)
    wvg_hnc, ssf_hnc, lfc_hnc, _ = run_stls_iet(rs, theta, "STLS-HNC")
    wvg_shnc, ssf_shnc, lfc_shnc, bf_shnc = run_stls_iet(rs, theta, "STLS-SHNC")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    axes[0].plot(wvg_stls, ssf_stls, label="STLS")
    axes[0].plot(wvg_hnc, ssf_hnc, "--", label="STLS-HNC")
    axes[0].plot(wvg_shnc, ssf_shnc, "-.", label="STLS-SHNC")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("S(x)")
    axes[0].set_title(rf"2D SSF, $r_s={rs:g}$, $\Theta={theta:g}$")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(wvg_stls, lfc_stls, label="STLS")
    axes[1].plot(wvg_hnc, lfc_hnc, "--", label="STLS-HNC")
    axes[1].plot(wvg_shnc, lfc_shnc, "-.", label="STLS-SHNC")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("G(x)")
    axes[1].set_title(rf"2D LFC, $r_s={rs:g}$, $\Theta={theta:g}$")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    axes[2].plot(wvg_shnc, bf_shnc, label=r"$B(x)/\beta U(x)$")
    axes[2].set_xlabel("x")
    axes[2].set_ylabel(r"$B(x)/\beta U(x)$")
    axes[2].set_title("Scaled-HNC prefactor")
    axes[2].grid(alpha=0.3)
    axes[2].legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
