import pytest
from numpy import linspace

from qupled.native import HF, Input


def test_hf_properties():
    scheme = HF(Input())
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "slfc")
    assert hasattr(scheme, "ssf")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "uint")
    assert excinfo.value.args[0] == "No data to compute the internal energy"


def test_hf_compute():
    inputs = Input()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "HF"
    inputs.chemical_potential = [-10, 10]
    inputs.matsubara = 128
    inputs.integral_error = 1.0e-5
    inputs.wave_vector_grid = linspace(0.0, 10, 100)
    inputs.threads = 1
    scheme = HF(inputs)
    scheme.compute()
    nx = inputs.wave_vector_grid.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.slfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(inputs.wave_vector_grid).size == nx
