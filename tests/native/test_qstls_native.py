import pytest
from numpy import linspace

from qupled.hf import DatabaseInfo
from qupled.native import Qstls, QstlsInput, Stls


@pytest.fixture
def database_info():
    dbInfo = DatabaseInfo()
    dbInfo.run_id = 1
    return dbInfo.to_native()


def test_qstls_properties():
    assert issubclass(Qstls, Stls)
    scheme = Qstls(QstlsInput())
    assert hasattr(scheme, "adr")


def test_qstls_compute(database_info):
    inputs = QstlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "QSTLS"
    inputs.chemical_potential = [-10, 10]
    inputs.matsubara = 32
    inputs.integral_error = 1.0e-5
    inputs.threads = 16
    inputs.wave_vector_grid = linspace(0.0, 5, 50)
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.database_info = database_info
    scheme = Qstls(inputs)
    scheme.compute()
    nx = inputs.wave_vector_grid.size
    assert nx >= 3
    assert scheme.adr.shape[0] == nx
    assert scheme.adr.shape[1] == inputs.matsubara
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.slfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(inputs.wave_vector_grid).size == nx
