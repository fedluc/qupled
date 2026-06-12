import json
from unittest import mock

import pytest

from qupled.mpi.worker_files import WorkerFiles


@pytest.fixture
def worker_files():
    worker_files = WorkerFiles()
    try:
        yield worker_files
    finally:
        worker_files.cleanup()


@pytest.mark.unit
def test_worker_files_owns_temporary_worker_directory():
    worker_files = WorkerFiles()
    worker_directory = worker_files.directory
    try:
        assert worker_directory.exists()
        assert worker_files.input_file == worker_directory / "input.json"
        assert worker_files.result_file == worker_directory / "results.json"
        assert worker_files.status_file == worker_directory / "status.json"
    finally:
        worker_files.cleanup()
    assert not worker_directory.exists()


@pytest.mark.unit
def test_worker_files_accepts_worker_directory(tmp_path):
    worker_files = WorkerFiles(tmp_path)
    assert worker_files.directory == tmp_path
    assert worker_files.input_file == tmp_path / "input.json"


@pytest.mark.unit
def test_write_inputs(worker_files):
    mock_inputs = mock.Mock()
    mock_inputs.to_dict.return_value = {"a": 1, "b": 2}
    worker_files.write_inputs(mock_inputs)
    with worker_files.input_file.open() as f:
        data = json.load(f)
    assert data == {"a": 1, "b": 2}
    mock_inputs.to_dict.assert_called_once()


@pytest.mark.unit
def test_read_inputs(worker_files):
    input_data = {"x": 10, "y": 20}
    with worker_files.input_file.open("w") as f:
        json.dump(input_data, f)
    mock_InputCls = mock.Mock()
    mock_InputCls.from_dict.return_value = "mocked_instance"
    result = worker_files.read_inputs(mock_InputCls)
    mock_InputCls.from_dict.assert_called_once_with(input_data)
    assert result == "mocked_instance"


@pytest.mark.unit
def test_write_results_writes_file_if_root(worker_files):
    mock_scheme = mock.Mock()
    mock_scheme.is_root = True
    mock_ResultCls = mock.Mock()
    mock_results_instance = mock.Mock()
    mock_ResultCls.return_value = mock_results_instance
    mock_results_instance.to_dict.return_value = {"foo": "bar"}
    worker_files.write_results(mock_scheme, mock_ResultCls)
    mock_ResultCls.assert_called_once_with()
    mock_results_instance.from_native.assert_called_once_with(mock_scheme)
    with worker_files.result_file.open() as f:
        data = json.load(f)
    assert data == {"foo": "bar"}
    mock_results_instance.to_dict.assert_called_once()


@pytest.mark.unit
def test_write_results_does_nothing_if_not_root(worker_files):
    mock_scheme = mock.Mock()
    mock_scheme.is_root = False
    mock_ResultCls = mock.Mock()
    worker_files.write_results(mock_scheme, mock_ResultCls)
    assert not worker_files.result_file.exists()
    mock_ResultCls.assert_not_called()


@pytest.mark.unit
def test_read_results(worker_files):
    input_data = {"x": 10, "y": 20}
    with worker_files.result_file.open("w") as f:
        json.dump(input_data, f)
    mock_ResultCls = mock.Mock()
    mock_ResultCls.from_dict.return_value = "mocked_instance"
    result = worker_files.read_results(mock_ResultCls)
    mock_ResultCls.from_dict.assert_called_once_with(input_data)
    assert result == "mocked_instance"


@pytest.mark.unit
def test_write_status_writes_file_if_root(worker_files):
    mock_scheme = mock.Mock()
    mock_scheme.is_root = True
    status_data = {"status": "done"}
    worker_files.write_status(mock_scheme, status_data)
    with worker_files.status_file.open() as f:
        data = json.load(f)
    assert data == status_data


@pytest.mark.unit
def test_write_status_does_nothing_if_not_root(worker_files):
    mock_scheme = mock.Mock()
    mock_scheme.is_root = False
    worker_files.write_status(mock_scheme, mock.ANY)
    assert not worker_files.status_file.exists()


@pytest.mark.unit
def test_read_status(worker_files):
    status_data = {"status": "ok", "code": 200}
    with worker_files.status_file.open("w") as f:
        json.dump(status_data, f)
    result = worker_files.read_status()
    assert result == status_data
