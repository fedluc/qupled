import os
from qupled.native import VSStls, VSStlsInput, Rpa


def test_vsstls_properties():
    assert issubclass(VSStls, Rpa)
    inputs = VSStlsInput()
    inputs.coupling = 1.0
    inputs.coupling_resolution = 0.1
    scheme = VSStls(inputs)
    assert hasattr(scheme, "free_energy_integrand")
    assert hasattr(scheme, "free_energy_grid")


def test_vsstls_compute():
    inputs = VSStlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "VSSTLS"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 5.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.output_frequency = 10
    inputs.coupling_resolution = 0.1
    inputs.degeneracy_resolution = 0.1
    inputs.error_alpha = 1.0e-3
    inputs.iterations_alpha = 50
    inputs.alpha = [0.5, 1.0]
    scheme = VSStls(inputs)
    scheme.compute()
    try:
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.ssf_HF.size == nx
        assert scheme.recovery == "recovery_rs1.000_theta1.000_VSSTLS.bin"
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recovery)
