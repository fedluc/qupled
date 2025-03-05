import pytest
from qupled.esa import ESA
from qupled.util import MPI


@pytest.fixture
def esa():
    return ESA()


@pytest.fixture
def esa_input():
    return ESA.Input(1.0, 1.0)


def test_default(esa):
    assert esa.hdf_file_name is None


def test_compute(esa, esa_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(ESA, ESA._compute.__name__)
    mock_save = mocker.patch.object(ESA, ESA._save.__name__)
    esa.compute(esa_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1
