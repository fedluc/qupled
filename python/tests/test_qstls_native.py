import os
import pytest
import glob
import set_path
import qupled.qupled as qp
import qupled.quantum as qpq

def test_qstls_properties():
    issubclass(qp.Qstls, qp.Stls)
    inputs = qpq.Qstls(1.0, 1.0).inputs
    scheme = qp.Qstls(inputs)
    assert hasattr(scheme, "adr")

def test_qstls_compute():
    inputs = qpq.Qstls(1.0, 1.0,
                       matsubara=32,
                       cutoff=5,
                       outputFrequency=2,
                       threads=16).inputs
    scheme = qp.Qstls(inputs)
    scheme.compute()
    try:
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.adr.shape[0] == nx
        assert scheme.adr.shape[1] == inputs.matsubara
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.ssfHF.size == nx
        assert scheme.recovery == "recovery_rs1.000_theta1.000_QSTLS.bin"
        assert os.path.isfile(scheme.recovery)
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        fixedFile = "adr_fixed_rs1.000_theta1.000_QSTLS.bin"
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recovery)
        if (os.path.isfile(fixedFile)) : os.remove(fixedFile)


def test_qstls_iet_properties():
    inputs = qpq.QstlsIet(1.0, 1.0, "QSTLS-HNC").inputs
    scheme = qp.Qstls(inputs)
    assert hasattr(scheme, "bf")

def test_qstls_iet_compute():
    ietSchemes = {"QSTLS-HNC",
                  "QSTLS-IOI",
                  "QSTLS-LCT"}
    for schemeName in ietSchemes:
        inputs = qpq.QstlsIet(10.0, 1.0, schemeName,
                              matsubara=16,
                              cutoff=5,
                              outputFrequency=2,
                              mixing=0.5,
                              threads=16).inputs
        scheme = qp.Qstls(inputs)
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
            if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recovery)
            fileNames = glob.glob("adr_fixed*.bin")
            for fileName in fileNames :
                os.remove(fileName)
