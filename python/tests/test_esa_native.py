import os
import pytest
import math
import set_path
import qupled.qupled as qp
import qupled.classic as qpc

def tolerance():
    return 1e-10

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



