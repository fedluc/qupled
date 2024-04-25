import os
import pytest
import set_path
import qupled.qupled as qp
import qupled.classic as qpc


def test_stls_properties():
    assert issubclass(qp.Stls, qp.Rpa)
    scheme = qp.Stls(qp.StlsInput())
    assert hasattr(scheme, "bf")


def test_stls_compute():
    inputs = qpc.Stls(1.0, 1.0).inputs
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
        inputs = qpc.StlsIet(
            10.0, 1.0, schemeName, matsubara=32, cutoff=5, outputFrequency=2, mixing=0.5
        ).inputs
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
