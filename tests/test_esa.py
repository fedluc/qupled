import pytest
from qupled.util import MPI
from qupled.classic import ESA


@pytest.fixture
def esa():
    return ESA()


@pytest.fixture
def esa_input():
    return ESA.Input(1.0, 1.0)


def test_default(esa):
    assert esa.hdfFileName == None


def test_compute(esa, esa_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(ESA, ESA._compute.__name__)
    mockSave = mocker.patch.object(ESA, ESA._save.__name__)
    esa.compute(esa_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1
