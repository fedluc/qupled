import sys
from unittest import mock

import pytest

from qupled.mpi import runtime
from qupled.mpi.worker_files import WorkerFiles


class MockSolver:
    pass


SOLVER_REFERENCE = f"{MockSolver.__module__}:{MockSolver.__qualname__}"


@pytest.fixture
def worker_files(tmp_path):
    return WorkerFiles(tmp_path)


@pytest.fixture
def mock_subprocess_run():
    with mock.patch("subprocess.run") as m:
        yield m


@pytest.fixture
def mock_shutil_which():
    with mock.patch("shutil.which") as m:
        yield m


@pytest.fixture
def mock_native():
    with mock.patch("qupled.mpi.runtime.native") as m:
        yield m


def expected_worker_command(worker_files):
    return [
        sys.executable,
        "-m",
        runtime._WORKER_MODULE,
        "--solver",
        SOLVER_REFERENCE,
        "--worker-directory",
        str(worker_files.directory),
    ]


@pytest.mark.unit
def test_launch_mpi_execution_calls_mpi(
    mock_subprocess_run, mock_shutil_which, mock_native, worker_files, mocker
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = True
    runtime._launch_mpi_execution(MockSolver, 4, worker_files)
    mock_subprocess_run.assert_called_once()
    command = mock_subprocess_run.call_args.args[0]
    assert command == [
        "mpiexec",
        "-n",
        "4",
        *expected_worker_command(worker_files),
    ]
    assert mock_subprocess_run.call_args.kwargs["check"] is True


@pytest.mark.unit
def test_launch_mpi_execution_serial_if_no_mpi(
    mock_subprocess_run, mock_shutil_which, mock_native, worker_files, capsys
):
    mock_shutil_which.return_value = None
    mock_native.uses_mpi = True
    runtime._launch_mpi_execution(MockSolver, 2, worker_files)
    mock_subprocess_run.assert_called_once()
    command = mock_subprocess_run.call_args.args[0]
    assert command == expected_worker_command(worker_files)
    assert mock_subprocess_run.call_args.kwargs["check"] is True
    captured = capsys.readouterr()
    assert "WARNING: Could not call MPI" in captured.out


@pytest.mark.unit
def test_launch_mpi_execution_serial_if_native_disabled(
    mock_subprocess_run,
    mock_shutil_which,
    mock_native,
    worker_files,
    capsys,
    mocker,
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = False
    runtime._launch_mpi_execution(MockSolver, 3, worker_files)
    mock_subprocess_run.assert_called_once()
    command = mock_subprocess_run.call_args.args[0]
    assert command == expected_worker_command(worker_files)
    assert mock_subprocess_run.call_args.kwargs["check"] is True
    captured = capsys.readouterr()
    assert "WARNING: Could not call MPI" in captured.out


@pytest.mark.unit
def test_launch_mpi_execution_raises_on_subprocess_error(
    mock_subprocess_run, mock_shutil_which, mock_native, worker_files, mocker
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = True
    mock_subprocess_run.side_effect = Exception("subprocess failed")
    with pytest.raises(Exception):
        runtime._launch_mpi_execution(MockSolver, 2, worker_files)


@pytest.mark.unit
def test_run_solver_owns_worker_files_and_returns_status_results(mocker):
    mock_inputs = mocker.Mock()
    mock_ResultCls = mocker.Mock()
    worker_files = mocker.Mock()
    worker_files.read_status.return_value = "mocked-status"
    worker_files.read_results.return_value = "mocked-results"
    worker_files_cls = mocker.patch(
        "qupled.mpi.runtime.WorkerFiles", return_value=worker_files
    )
    launch_mpi_execution = mocker.patch("qupled.mpi.runtime._launch_mpi_execution")
    workflow = mocker.Mock()
    workflow.attach_mock(worker_files.write_inputs, "write_inputs")
    workflow.attach_mock(launch_mpi_execution, "launch_mpi_execution")
    status, results = runtime.run_solver(MockSolver, mock_inputs, 4, mock_ResultCls)
    worker_files_cls.assert_called_once_with()
    worker_files.write_inputs.assert_called_once_with(mock_inputs)
    launch_mpi_execution.assert_called_once_with(MockSolver, 4, worker_files)
    assert workflow.mock_calls == [
        mock.call.write_inputs(mock_inputs),
        mock.call.launch_mpi_execution(MockSolver, 4, worker_files),
    ]
    worker_files.read_status.assert_called_once_with()
    worker_files.read_results.assert_called_once_with(mock_ResultCls)
    worker_files.cleanup.assert_called_once_with()
    assert status == "mocked-status"
    assert results == "mocked-results"


@pytest.mark.unit
def test_run_solver_cleans_up_on_failure(mocker):
    worker_files = mocker.Mock()
    mocker.patch("qupled.mpi.runtime.WorkerFiles", return_value=worker_files)
    launch_mpi_execution = mocker.patch("qupled.mpi.runtime._launch_mpi_execution")
    launch_mpi_execution.side_effect = RuntimeError
    with pytest.raises(RuntimeError):
        runtime.run_solver(MockSolver, mocker.Mock(), 2, mocker.Mock())
    worker_files.cleanup.assert_called_once_with()
