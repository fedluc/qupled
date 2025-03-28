import pytest
from unittest import mock
from qupled.util import DataBase, MPI


@pytest.fixture
def db_handler():
    with mock.patch("qupled.util.database.DataBaseHandler") as MockDBHandler:
        yield MockDBHandler


def test_inspect_runs(db_handler):
    db_handler.return_value.inspect_runs.return_value = {"run1": "data1"}
    result = DataBase.inspect_runs("test_db")
    assert result == {"run1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_run(db_handler):
    db_handler.return_value.get_run.return_value = {
        "input1": "data1",
        "result1": "data2",
    }
    result = DataBase.read_run(1, "test_db", ["input1"], ["result1"])
    assert result == {"input1": "data1", "result1": "data2"}
    db_handler.assert_called_once_with("test_db")


def test_read_inputs(db_handler):
    db_handler.return_value.get_inputs.return_value = {"input1": "data1"}
    result = DataBase.read_inputs(1, "test_db", ["input1"])
    assert result == {"input1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_results(db_handler):
    db_handler.return_value.get_results.return_value = {"result1": "data1"}
    result = DataBase.read_results(1, "test_db", ["result1"])
    assert result == {"result1": "data1"}
    db_handler.assert_called_once_with("test_db")


@pytest.fixture
def mpi_native():
    with mock.patch("qupled.util.native.MPI") as MockMPI:
        yield MockMPI


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


def test_record_time_decorator(mpi_native):
    mpi_native.return_value.timer.side_effect = [0, 3600]

    @MPI.record_time
    def test_func():
        pass

    with mock.patch("builtins.print") as mock_print:
        test_func()
        mock_print.assert_called_once_with("Elapsed time: 1 h, 0 m, 0 s.")
