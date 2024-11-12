import os
import pytest
import qupled.qupled as qp


def test_stls_properties():
    assert issubclass(qp.Stls, qp.Rpa)
    scheme = qp.Stls(qp.StlsInput())
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "slfc")
    assert hasattr(scheme, "ssf")
    assert hasattr(scheme, "ssfHF")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "uInt")
    assert excinfo.value.args[0] == "No data to compute the internal energy"
    assert hasattr(scheme, "wvg")
    assert hasattr(scheme, "recovery")
    assert hasattr(scheme, "bf")


def test_stls_compute():
    inputs = qp.StlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "STLS"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.outputFrequency = 10
    scheme = qp.Stls(inputs)
    scheme.compute()
    try:
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.ssfHF.size == nx
        assert scheme.recovery == "recovery_rs1.000_theta1.000_STLS.bin"
        assert os.path.isfile(scheme.recovery)
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recovery)


def test_stls_iet_compute():
    ietSchemes = {"STLS-HNC", "STLS-IOI", "STLS-LCT"}
    for schemeName in ietSchemes:
        inputs = qp.StlsInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = schemeName
        inputs.chemicalPotential = [-10, 10]
        inputs.cutoff = 5.0
        inputs.matsubara = 32
        inputs.resolution = 0.1
        inputs.intError = 1.0e-5
        inputs.threads = 1
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        inputs.outputFrequency = 2
        scheme = qp.Stls(inputs)
        scheme.compute()
        try:
            nx = scheme.wvg.size
            assert nx >= 3
            assert scheme.idr.shape[0] == nx
            assert scheme.idr.shape[1] == inputs.matsubara
            assert scheme.sdr.size == nx
            assert scheme.slfc.size == nx
            assert scheme.ssf.size == nx
            assert scheme.ssfHF.size == nx
            recovery = "recovery_rs10.000_theta1.000_" + schemeName + ".bin"
            assert scheme.recovery == recovery
            assert os.path.isfile(scheme.recovery)
            assert scheme.rdf(scheme.wvg).size == nx
        finally:
            if os.path.isfile(scheme.recovery):
                os.remove(scheme.recovery)
