import os
import pytest
import set_path
import qupled.qupled as qp
import qupled.classic as qpc

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



