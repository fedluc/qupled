import pytest

from qupled.native import Stls, Rpa, StlsInput


def test_stls_properties():
    assert issubclass(Stls, Rpa)


def test_stls_compute():
    inputs = StlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "STLS"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    scheme = Stls(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.lfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(scheme.wvg).size == nx
