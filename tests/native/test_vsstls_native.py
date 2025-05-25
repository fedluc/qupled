import pytest
from qupled.native import VSStls, VSStlsInput, Rpa


def test_vsstls_properties():
    inputs = VSStlsInput()
    inputs.coupling = 1.0
    inputs.coupling_resolution = 0.1
    scheme = VSStls(inputs)
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "lfc")
    assert hasattr(scheme, "ssf")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "uint")
    assert excinfo.value.args[0] == "No data to compute the internal energy"
    assert hasattr(scheme, "wvg")
    assert hasattr(scheme, "error")
    assert hasattr(scheme, "free_energy_integrand")
    assert hasattr(scheme, "free_energy_grid")


def test_vsstls_compute():
    inputs = VSStlsInput()
    inputs.coupling = 0.1
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
    inputs.coupling_resolution = 0.1
    inputs.degeneracy_resolution = 0.1
    inputs.error_alpha = 1.0e-3
    inputs.iterations_alpha = 50
    inputs.alpha = [0.5, 1.0]
    scheme = VSStls(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.lfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(scheme.wvg).size == nx
