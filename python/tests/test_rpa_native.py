import os
import pytest
import set_path
import qupled.qupled as qp
import qupled.classic as qpc

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
    assert round(scheme.uInt, 5) == -0.52298



