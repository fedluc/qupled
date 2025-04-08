from unittest.mock import PropertyMock
import pytest
import numpy as np
from qupled.database import DataBaseHandler
import qupled.hf as hf
import qupled.native as native


@pytest.fixture
def inputs():
    return hf.Input(coupling=1.0, degeneracy=2.0)


@pytest.fixture
def results():
    return hf.Result()


@pytest.fixture
def scheme(mocker):
    scheme = hf.HF()
    scheme.db_handler = mocker.Mock()
    return scheme


def test_hf_initialization():
    scheme = hf.HF()
    assert scheme.inputs is None
    assert isinstance(scheme.results, hf.Result)
    assert isinstance(scheme.db_handler, DataBaseHandler)
    assert scheme.native_scheme_cls == native.HF
    assert isinstance(scheme.native_inputs, native.Input)
    assert scheme.native_scheme_status is None


def test_run_id(scheme):
    run_id = "run_id"
    scheme.db_handler.run_id = run_id
    assert scheme.run_id == run_id


def test_compute(scheme, inputs, mocker):
    add_run_to_database = mocker.patch.object(scheme, "_add_run_to_database")
    compute_native = mocker.patch.object(scheme, "_compute_native")
    save = mocker.patch.object(scheme, "_save")
    scheme.compute(inputs)
    assert scheme.inputs is not None
    add_run_to_database.assert_called_once()
    compute_native.assert_called_once()
    save.assert_called_once()


def test_add_run_to_database(scheme, mocker):
    mocker.patch.object(hf.HF, "run_id", new_callable=PropertyMock).return_value = (
        "mocked-run-id"
    )
    scheme.inputs = mocker.Mock()
    scheme._add_run_to_database()
    scheme.db_handler.insert_run.assert_called_once_with(scheme.inputs)
    assert scheme.inputs.database_info.run_id == scheme.run_id


def test_compute_native(scheme, inputs, mocker):
    to_native = mocker.patch("qupled.hf.Input.to_native")
    native_scheme = mocker.Mock()
    native_scheme_cls = mocker.patch.object(
        scheme, "native_scheme_cls", return_value=native_scheme
    )
    from_native = mocker.patch("qupled.hf.Result.from_native")
    native_scheme.compute.return_value = "mocked-status"
    scheme.inputs = inputs
    scheme._compute_native()
    to_native.assert_called_once_with(scheme.native_inputs)
    native_scheme_cls.assert_called_once_with(scheme.native_inputs)
    native_scheme.compute.assert_called_once()
    assert scheme.native_scheme_status == "mocked-status"
    from_native.assert_called_once_with(native_scheme)


def test_save(scheme, results, mocker):
    scheme.results = results
    scheme.native_scheme_status = mocker.Mock()
    scheme._save()
    scheme.db_handler.update_run_status.assert_called_once_with(
        scheme.native_scheme_status
    )
    scheme.db_handler.insert_results.assert_called_once_with(scheme.results.__dict__)


def test_compute_rdf_with_default_grid(scheme, results, mocker):
    compute_rdf = mocker.patch("qupled.hf.Result.compute_rdf")
    scheme.results = results
    scheme.compute_rdf()
    compute_rdf.assert_called_once_with(None)
    scheme.db_handler.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        }
    )


def test_compute_rdf_with_custom_grid(scheme, results, mocker):
    compute_rdf = mocker.patch("qupled.hf.Result.compute_rdf")
    scheme.results = results
    rdf_grid = np.array([1, 2, 3])
    scheme.compute_rdf(rdf_grid)
    compute_rdf.assert_called_once_with(rdf_grid)
    scheme.db_handler.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        }
    )


def test_compute_rdf_without_results(scheme):
    scheme.results = None
    scheme.compute_rdf()
    scheme.db_handler.insert_results.assert_not_called()


def test_input_initialization(inputs):
    assert inputs.coupling == 1.0
    assert inputs.degeneracy == 2.0
    assert inputs.chemical_potential == [-10.0, 10.0]
    assert inputs.cutoff == 10.0
    assert inputs.frequency_cutoff == 10.0
    assert inputs.integral_error == 1.0e-5
    assert inputs.integral_strategy == "full"
    assert inputs.matsubara == 128
    assert inputs.resolution == 0.1
    assert inputs.threads == 1
    assert inputs.theory == "HF"


def test_input_to_native(mocker, inputs):
    native_input = mocker.Mock()
    inputs.to_native(native_input)
    assert native_input.coupling == 1.0
    assert native_input.degeneracy == 2.0


def test_result_initialization(results):
    assert results.idr is None
    assert results.rdf is None
    assert results.rdf_grid is None
    assert results.sdr is None
    assert results.slfc is None
    assert results.ssf is None
    assert results.uint is None
    assert results.wvg is None


def test_result_from_native(mocker, results):
    native_scheme = mocker.Mock()
    native_scheme.idr = np.array([1, 2, 3])
    results.from_native(native_scheme)
    assert np.array_equal(results.idr, np.array([1, 2, 3]))


def test_result_compute_rdf_with_default_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    results.wvg = np.array([1, 2, 3])
    results.ssf = np.array([4, 5, 6])
    native_compute_rdf.return_value = np.array([7, 8, 9])
    results.compute_rdf()
    assert results.rdf is not None
    assert results.rdf_grid is not None
    native_compute_rdf.assert_called_once()
    assert np.array_equal(results.rdf, np.array([7, 8, 9]))
    assert np.array_equal(results.rdf_grid, np.arange(0.0, 10.0, 0.01))


def test_result_compute_rdf_with_custom_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    rdf_grid = np.array([0, 1, 2])
    results.wvg = np.array([1, 2, 3])
    results.ssf = np.array([4, 5, 6])
    native_compute_rdf.return_value = np.array([7, 8, 9])
    results.compute_rdf(rdf_grid)
    assert results.rdf is not None
    assert results.rdf_grid is not None
    native_compute_rdf.assert_called_once()
    assert np.array_equal(results.rdf, np.array([7, 8, 9]))
    assert np.array_equal(results.rdf_grid, rdf_grid)
