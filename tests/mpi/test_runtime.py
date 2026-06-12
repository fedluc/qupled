import json
import sys
from unittest import mock

import pytest

from qupled.mpi import runtime


class LauncherSolver:
    pass


SOLVER_REFERENCE = f"{LauncherSolver.__module__}:{LauncherSolver.__qualname__}"


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


@pytest.mark.unit
def test_launch_mpi_execution_calls_mpi(
    mock_subprocess_run, mock_shutil_which, mock_native, mocker
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = True
    runtime.launch_mpi_execution(LauncherSolver, 4)
    mock_subprocess_run.assert_called_once_with(
        [
            "mpiexec",
            "-n",
            "4",
            sys.executable,
            "-m",
            runtime._WORKER_MODULE,
            "--solver",
            SOLVER_REFERENCE,
        ],
        check=True,
    )


@pytest.mark.unit
def test_launch_mpi_execution_serial_if_no_mpi(
    mock_subprocess_run, mock_shutil_which, mock_native, capsys
):
    mock_shutil_which.return_value = None
    mock_native.uses_mpi = True
    runtime.launch_mpi_execution(LauncherSolver, 2)
    mock_subprocess_run.assert_called_once_with(
        [
            sys.executable,
            "-m",
            runtime._WORKER_MODULE,
            "--solver",
            SOLVER_REFERENCE,
        ],
        check=True,
    )
    captured = capsys.readouterr()
    assert "WARNING: Could not call MPI" in captured.out


@pytest.mark.unit
def test_launch_mpi_execution_serial_if_native_disabled(
    mock_subprocess_run, mock_shutil_which, mock_native, capsys, mocker
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = False
    runtime.launch_mpi_execution(LauncherSolver, 3)
    mock_subprocess_run.assert_called_once_with(
        [
            sys.executable,
            "-m",
            runtime._WORKER_MODULE,
            "--solver",
            SOLVER_REFERENCE,
        ],
        check=True,
    )
    captured = capsys.readouterr()
    assert "WARNING: Could not call MPI" in captured.out


@pytest.mark.unit
def test_launch_mpi_execution_raises_on_subprocess_error(
    mock_subprocess_run, mock_shutil_which, mock_native, mocker
):
    mock_shutil_which.return_value = mocker.ANY
    mock_native.uses_mpi = True
    mock_subprocess_run.side_effect = Exception("subprocess failed")
    with pytest.raises(Exception):
        runtime.launch_mpi_execution(LauncherSolver, 2)


@pytest.mark.unit
def test_write_inputs(tmp_path, mocker):
    test_file = tmp_path / "input.json"
    mocker.patch("qupled.mpi.runtime.INPUT_FILE", test_file)
    mock_inputs = mock.Mock()
    mock_inputs.to_dict.return_value = {"a": 1, "b": 2}
    runtime.write_inputs(mock_inputs)
    with test_file.open() as f:
        data = json.load(f)
    assert data == {"a": 1, "b": 2}
    mock_inputs.to_dict.assert_called_once()


@pytest.mark.unit
def test_read_inputs(tmp_path, mocker):
    test_file = tmp_path / "input.json"
    mocker.patch("qupled.mpi.runtime.INPUT_FILE", test_file)
    input_data = {"x": 10, "y": 20}
    with test_file.open("w") as f:
        json.dump(input_data, f)
    mock_InputCls = mock.Mock()
    mock_InputCls.from_dict.return_value = "mocked_instance"
    result = runtime.read_inputs(mock_InputCls)
    mock_InputCls.from_dict.assert_called_once_with(input_data)
    assert result == "mocked_instance"


@pytest.mark.unit
def test_write_results_writes_file_if_root(tmp_path, mocker):
    test_file = tmp_path / "results.json"
    mocker.patch("qupled.mpi.runtime.RESULT_FILE", test_file)
    mock_scheme = mock.Mock()
    mock_scheme.is_root = True
    mock_ResultCls = mock.Mock()
    mock_results_instance = mock.Mock()
    mock_ResultCls.return_value = mock_results_instance
    mock_results_instance.to_dict.return_value = {"foo": "bar"}
    runtime.write_results(mock_scheme, mock_ResultCls)
    mock_ResultCls.assert_called_once_with()
    mock_results_instance.from_native.assert_called_once_with(mock_scheme)
    with test_file.open() as f:
        data = json.load(f)
    assert data == {"foo": "bar"}
    mock_results_instance.to_dict.assert_called_once()


@pytest.mark.unit
def test_write_results_does_nothing_if_not_root(tmp_path, mocker):
    test_file = tmp_path / "results.json"
    mocker.patch("qupled.mpi.runtime.RESULT_FILE", test_file)
    mock_scheme = mock.Mock()
    mock_scheme.is_root = False
    mock_ResultCls = mock.Mock()
    runtime.write_results(mock_scheme, mock_ResultCls)
    assert not test_file.exists()
    mock_ResultCls.assert_not_called()


@pytest.mark.unit
def test_read_results(tmp_path, mocker):
    test_file = tmp_path / "input.json"
    mocker.patch("qupled.mpi.runtime.RESULT_FILE", test_file)
    input_data = {"x": 10, "y": 20}
    with test_file.open("w") as f:
        json.dump(input_data, f)
    mock_InputCls = mock.Mock()
    mock_InputCls.from_dict.return_value = "mocked_instance"
    result = runtime.read_results(mock_InputCls)
    mock_InputCls.from_dict.assert_called_once_with(input_data)
    assert result == "mocked_instance"


@pytest.mark.unit
def test_write_status_writes_file_if_root(tmp_path, mocker):
    test_file = tmp_path / "status.json"
    mocker.patch("qupled.mpi.runtime.STATUS_FILE", test_file)
    mock_scheme = mock.Mock()
    mock_scheme.is_root = True
    status_data = {"status": "done"}
    runtime.write_status(mock_scheme, status_data)
    with test_file.open() as f:
        data = json.load(f)
    assert data == status_data


@pytest.mark.unit
def test_write_status_does_nothing_if_not_root(tmp_path, mocker):
    test_file = tmp_path / "status.json"
    mocker.patch("qupled.mpi.runtime.STATUS_FILE", test_file)
    mock_scheme = mock.Mock()
    mock_scheme.is_root = False
    runtime.write_status(mock_scheme, mocker.ANY)
    assert not test_file.exists()


@pytest.mark.unit
def test_read_status(tmp_path, mocker):
    test_file = tmp_path / "status.json"
    mocker.patch("qupled.mpi.runtime.STATUS_FILE", test_file)
    status_data = {"status": "ok", "code": 200}
    with test_file.open("w") as f:
        json.dump(status_data, f)
    result = runtime.read_status()
    assert result == status_data


@pytest.mark.unit
def test_clean_files_removes_existing_files(tmp_path, mocker):
    input_file = tmp_path / "input.json"
    result_file = tmp_path / "results.json"
    status_file = tmp_path / "status.json"
    for f in [input_file, result_file, status_file]:
        f.write_text("test")
    mocker.patch("qupled.mpi.runtime.INPUT_FILE", input_file)
    mocker.patch("qupled.mpi.runtime.RESULT_FILE", result_file)
    mocker.patch("qupled.mpi.runtime.STATUS_FILE", status_file)
    runtime.clean_files()
    assert not input_file.exists()
    assert not result_file.exists()
    assert not status_file.exists()


@pytest.mark.unit
def test_clean_files_does_nothing_if_files_do_not_exist(tmp_path, mocker):
    input_file = tmp_path / "input.json"
    result_file = tmp_path / "results.json"
    status_file = tmp_path / "status.json"
    for f in [input_file, result_file, status_file]:
        if f.exists():
            f.unlink()
    mocker.patch("qupled.mpi.runtime.INPUT_FILE", input_file)
    mocker.patch("qupled.mpi.runtime.RESULT_FILE", result_file)
    mocker.patch("qupled.mpi.runtime.STATUS_FILE", status_file)
    runtime.clean_files()
    assert not input_file.exists()
    assert not result_file.exists()
    assert not status_file.exists()
