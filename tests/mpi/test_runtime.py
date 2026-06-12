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
    runtime.launch_mpi_execution(MockSolver, 4, worker_files)
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
    runtime.launch_mpi_execution(MockSolver, 2, worker_files)
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
    runtime.launch_mpi_execution(MockSolver, 3, worker_files)
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
        runtime.launch_mpi_execution(MockSolver, 2, worker_files)
