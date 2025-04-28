import pytest
from numpy import linspace

from qupled.hf import DatabaseInfo
from qupled.native import Qstls, QstlsIet, QstlsIetInput


@pytest.fixture
def database_info():
    dbInfo = DatabaseInfo()
    dbInfo.run_id = 1
    return dbInfo.to_native()


def test_qstls_properties():
    assert issubclass(QstlsIet, Qstls)
    scheme = QstlsIet(QstlsIetInput())
    assert hasattr(scheme, "bf")


def test_qstls_iet_compute(database_info):
    iet_schemes = {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = QstlsIetInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = scheme_name
        inputs.chemical_potential = [-10, 10]
        inputs.matsubara = 16
        inputs.integral_error = 1.0e-5
        inputs.threads = 16
        inputs.wave_vector_grid = linspace(0.0, 5, 50)
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        inputs.database_info = database_info
        scheme = QstlsIet(inputs)
        scheme.compute()
        nx = inputs.wave_vector_grid.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.bf.size == nx
        assert scheme.rdf(inputs.wave_vector_grid).size == nx
