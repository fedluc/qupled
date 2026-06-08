import matplotlib.pyplot as plt
import numpy as np

from qupled.schemes import esa, recstls, stls

RS = 10.0
THETA = 1.0
CUTOFF = 10.0
RESOLUTION = 0.1
RDF_MAX = 10.0
RDF_STEP = 0.01


def main():
    print("######### Solving the ESA scheme #########")
    esa_scheme = esa.Solver()
    esa_inputs = esa.Input(
        coupling=RS,
        degeneracy=THETA,
        cutoff=CUTOFF,
        resolution=RESOLUTION,
    )
    esa_scheme.compute(esa_inputs)
    rdf_grid = np.arange(0.0, RDF_MAX + 0.5 * RDF_STEP, RDF_STEP)
    esa_scheme.results.compute_rdf(esa_scheme.inputs.dimension.value, rdf_grid)

    print("######### Solving the STLS scheme #########")
    stls_scheme = stls.Solver()
    stls_inputs = stls.Input(
        coupling=RS,
        degeneracy=THETA,
        cutoff=CUTOFF,
        resolution=RESOLUTION,
        mixing=0.2,
        iterations=200,
        error=1.0e-5,
    )
    stls_scheme.compute(stls_inputs)
    stls_scheme.results.compute_rdf(stls_scheme.inputs.dimension.value, rdf_grid)

    print("######### Solving the REC-STLS scheme #########")
    rec_scheme = recstls.Solver()
    rec_inputs = recstls.Input(
        coupling=RS,
        degeneracy=THETA,
        cutoff=CUTOFF,
        resolution=RESOLUTION,
        rdf_grid=esa_scheme.results.rdf_grid,
        rdf=esa_scheme.results.rdf,
    )
    rec_scheme.compute(rec_inputs)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(rec_scheme.results.wvg, rec_scheme.results.ssf_input, label=r"$S_{\rm rec}$")
    ax1.plot(rec_scheme.results.wvg, rec_scheme.results.ssf, label=r"$S_{\rm out}$")
    ax1.plot(stls_scheme.results.wvg, stls_scheme.results.ssf, "-.", label="STLS")
    ax1.plot(esa_scheme.results.wvg, esa_scheme.results.ssf, "--", label="ESA")
    ax1.set_xlabel(r"$q / k_F$")
    ax1.set_ylabel(r"$S(q)$")
    ax1.set_title(fr"Structure factors, $r_s={RS:g}$, $\theta={THETA:g}$")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(rec_inputs.rdf_grid, rec_inputs.rdf, label=r"$g_{\rm rec}$ input")
    ax2.plot(stls_scheme.results.rdf_grid, stls_scheme.results.rdf, "-.", label="STLS RDF")
    ax2.plot(esa_scheme.results.rdf_grid, esa_scheme.results.rdf, "--", label="ESA RDF")
    ax2.set_xlabel(r"$r k_F$")
    ax2.set_ylabel(r"$g(r)$")
    ax2.set_title("Input reconstructed RDF")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
