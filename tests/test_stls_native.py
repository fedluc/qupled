import os
import pytest
from qupled.native import Stls, Rpa, StlsInput


def test_stls_properties():
    assert issubclass(Stls, Rpa)
    scheme = Stls(StlsInput())
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "slfc")
    assert hasattr(scheme, "ssf")
    assert hasattr(scheme, "ssf_HF")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "internal_energy")
    assert excinfo.value.args[0] == "No data to compute the internal energy"
    assert hasattr(scheme, "wvg")
    assert hasattr(scheme, "recovery")
    assert hasattr(scheme, "bf")


def test_stls_compute():
    inputs = StlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "STLS"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.output_frequency = 10
    scheme = Stls(inputs)
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
        assert scheme.recovery == "recovery_rs1.000_theta1.000_STLS.bin"
        assert os.path.isfile(scheme.recovery)
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recovery)


def test_stls_iet_compute():
    iet_schemes = {"STLS-HNC", "STLS-IOI", "STLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = StlsInput()
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
        inputs.output_frequency = 2
        scheme = Stls(inputs)
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
            recovery = "recovery_rs10.000_theta1.000_" + scheme_name + ".bin"
            assert scheme.recovery == recovery
            assert os.path.isfile(scheme.recovery)
            assert scheme.rdf(scheme.wvg).size == nx
        finally:
            if os.path.isfile(scheme.recovery):
                os.remove(scheme.recovery)
