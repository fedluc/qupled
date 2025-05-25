import pytest
from qupled.native import Stls, StlsIet, StlsIetInput


def test_stls_properties():
    scheme = StlsIet(StlsIetInput())
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "lfc")
    assert hasattr(scheme, "ssf")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "uint")
    assert excinfo.value.args[0] == "No data to compute the internal energy"
    assert hasattr(scheme, "wvg")
    assert hasattr(scheme, "error")
    assert hasattr(scheme, "bf")


def test_stls_iet_compute():
    iet_schemes = {"STLS-HNC", "STLS-IOI", "STLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = StlsIetInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = scheme_name
        inputs.chemical_potential = [-10, 10]
        inputs.cutoff = 5.0
        inputs.matsubara = 32
        inputs.resolution = 0.1
        inputs.integral_error = 1.0e-5
        inputs.threads = 1
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        scheme = StlsIet(inputs)
        scheme.compute()
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.lfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.bf.size == nx
        assert scheme.rdf(scheme.wvg).size == nx
