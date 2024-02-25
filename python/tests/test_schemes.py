import os
import pytest
import math
import numpy as np
import set_path
import qupled.qupled as qp
import qupled.classic as qpc
import qupled.quantum as qpq

def tolerance():
    return 1e-10

def test_rpa_properties():
    scheme = qp.Rpa(qp.RpaInput())
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

def test_rpa_compute():
    inputs = qpc.Rpa(1.0, 1.0).inputs
    scheme = qp.Rpa(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.slfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.ssfHF.size == nx
    assert scheme.recovery == ""
    assert scheme.rdf(scheme.wvg).size == nx
    assert math.isclose(scheme.uInt, -0.5229843450, rel_tol=tolerance())

def test_esa_properties():
    issubclass(qp.ESA, qp.Rpa)

def test_esa_compute():
    inputs = qpc.ESA(1.0, 1.0).inputs
    scheme = qp.ESA(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.slfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.ssfHF.size == nx
    assert scheme.recovery == ""
    assert scheme.rdf(scheme.wvg).size == nx
    assert math.isclose(scheme.uInt, -0.4801292597, rel_tol=tolerance())
                                     
def test_stls_properties():
    issubclass(qp.Stls, qp.Rpa)
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
        assert math.isclose(scheme.uInt, -0.4862983323, rel_tol=tolerance())
    finally:
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recovery)

def test_stls_iet_compute():
    ietSchemes = {"STLS-HNC" : -0.07132783536,
                  "STLS-IOI" : -0.07132737479,
                  "STLS-LCT" : -0.07172490005}
    for schemeName, uInt in ietSchemes.items():
        inputs = qpc.StlsIet(10.0, 1.0, schemeName,
                             matsubara=32,
                             cutoff=5,
                             outputFrequency=2,
                             mixing=0.5).inputs
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
            assert math.isclose(scheme.uInt, uInt, rel_tol=tolerance())
        finally:
            if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recovery)


def test_vsstls_properties():
    issubclass(qp.VSStls, qp.Rpa)
    inputs = qpc.VSStls(1.0, 1.0).inputs
    scheme = qp.VSStls(inputs)
    assert hasattr(scheme, "freeEnergyIntegrand")
    assert hasattr(scheme, "freeEnergyGrid")
    
def test_stls_compute():
    inputs = qpc.VSStls(1.0, 1.0,
                        couplingResolution=0.1,
                        degeneracyResolution=0.1,
                        cutoff=5).inputs
    scheme = qp.VSStls(inputs)
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
        assert scheme.recovery == ""
        assert scheme.rdf(scheme.wvg).size == nx
        assert math.isclose(scheme.uInt, -0.5223739856, rel_tol=tolerance())
    finally:
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recover)


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
        assert math.isclose(scheme.uInt, -0.4873131719, rel_tol=tolerance())
    finally:
        fixedFile = "adr_fixed_rs1.000_theta1.000_QSTLS.bin"
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recovery)
        if (os.path.isfile(fixedFile)) : os.remove(fixedFile)



