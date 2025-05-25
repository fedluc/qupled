import pytest

from qupled.hf import DatabaseInfo
from qupled.native import Qstls, QstlsInput, Stls


@pytest.fixture
def database_info():
    dbInfo = DatabaseInfo()
    dbInfo.run_id = 1
    return dbInfo.to_native()


def test_qstls_properties():
    scheme = Qstls(QstlsInput())
    assert hasattr(scheme, "idr")
    assert hasattr(scheme, "sdr")
    assert hasattr(scheme, "lfc")
    assert hasattr(scheme, "ssf")
    with pytest.raises(RuntimeError) as excinfo:
        hasattr(scheme, "uint")
    assert excinfo.value.args[0] == "No data to compute the internal energy"
    assert hasattr(scheme, "wvg")
    assert hasattr(scheme, "error")


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
    assert scheme.lfc.shape[0] == nx
    assert scheme.lfc.shape[1] == inputs.matsubara
    assert scheme.idr.shape[0] == nx
    assert scheme.idr.shape[1] == inputs.matsubara
    assert scheme.sdr.size == nx
    assert scheme.ssf.size == nx
    assert scheme.rdf(scheme.wvg).size == nx
