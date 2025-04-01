import pytest
from qupled.mpi import MPI


@pytest.fixture
def mpi_native(mocker):
    yield mocker.patch("qupled.native.MPI")


def test_mpi_rank(mpi_native):
    mpi_native.return_value.rank.return_value = 0
    mpi = MPI()
    assert mpi.rank() == 0


def test_mpi_is_root(mpi_native):
    mpi_native.return_value.is_root.return_value = True
    mpi = MPI()
    assert mpi.is_root() is True


def test_mpi_barrier(mpi_native):
    mpi = MPI()
    mpi.barrier()
    mpi_native.return_value.barrier.assert_called_once()


def test_mpi_timer(mpi_native):
    mpi_native.return_value.timer.return_value = 123.456
    mpi = MPI()
    assert mpi.timer() == 123.456


def test_run_only_on_root_decorator(mpi_native):
    mpi_native.return_value.is_root.return_value = True

    @MPI.run_only_on_root
    def test_func():
        return "Executed"

    assert test_func() == "Executed"


def test_synchronize_ranks_decorator(mpi_native):
    @MPI.synchronize_ranks
    def test_func():
        pass

    test_func()
    mpi_native.return_value.barrier.assert_called_once()


def test_record_time_decorator(mpi_native, mocker):
    mpi_native.return_value.timer.side_effect = [0, 3600]

    @MPI.record_time
    def test_func():
        pass

    mock_print = mocker.patch("builtins.print")
    test_func()
    mock_print.assert_called_once_with("Elapsed time: 1 h, 0 m, 0 s.")
