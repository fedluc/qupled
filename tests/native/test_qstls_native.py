import pytest

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
    inputs.cutoff = 5.0
    inputs.matsubara = 32
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 16
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.database_info = database_info
    scheme = Qstls(inputs)
    scheme.compute()
    nx = scheme.wvg.size
    assert nx >= 3
    assert scheme.adr.shape[0] == nx
    assert scheme.adr.shape[1] == inputs.matsubara
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.slfc.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(scheme.wvg).size == nx


def test_qstls_iet_compute(database_info):
    iet_schemes = {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = QstlsInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = scheme_name
        inputs.chemical_potential = [-10, 10]
        inputs.cutoff = 5.0
        inputs.matsubara = 16
        inputs.resolution = 0.1
        inputs.integral_error = 1.0e-5
        inputs.threads = 16
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        inputs.database_info = database_info
        scheme = Qstls(inputs)
        scheme.compute()
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.rdf(scheme.wvg).size == nx
