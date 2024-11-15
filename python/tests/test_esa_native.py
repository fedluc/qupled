import os
import pytest
from qupled import native


def test_esa_properties():
    assert issubclass(native.ESA, native.Rpa)
    scheme = native.ESA(native.RpaInput())
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


def test_esa_compute():
    inputs = native.RpaInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "RPA"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    scheme = native.ESA(inputs)
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
