import pytest

from qupled.native import HF, Rpa, Input


def test_rpa_properties():
    assert issubclass(Rpa, HF)


def test_rpa_compute():
    inputs = Input()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "RPA"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 1
    scheme = Rpa(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.lfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(scheme.wvg).size == nx
