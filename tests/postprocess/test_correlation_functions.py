import pytest
import numpy as np

from qupled.postprocess.correlation_functions import compute_rdf, compute_itcf


@pytest.fixture
def db_handler(mocker):
    return mocker.patch("qupled.postprocess.correlation_functions.DataBaseHandler")


@pytest.fixture
def scheme_tables(db_handler):
    return db_handler.return_value.scheme_tables


@pytest.fixture
def read_run(mocker):
    return mocker.patch("qupled.postprocess.correlation_functions.DataBase.read_run")


@pytest.mark.unit
def test_compute_rdf_with_default_grid(mocker, read_run, db_handler, scheme_tables):
    mock_data = mocker.Mock()
    mock_data.inputs = {"dimension": {"value": "_3D"}}
    mock_data.results = {"wvg": np.array([1.0, 2.0]), "ssf": np.array([0.5, 0.6])}
    read_run.return_value = mock_data
    mock_result = mocker.Mock()
    mock_result.rdf = np.array([1.0, 2.0])
    mock_result.rdf_grid = np.array([0.0, 0.1])
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_result)
    mock_dimension = mocker.Mock()
    mocker.patch(
        "qupled.util.dimension.Dimension.from_dict", return_value=mock_dimension
    )
    compute_rdf(run_id=1)
    mock_result.compute_rdf.assert_called_once_with(mock_dimension, None)
    scheme_tables.insert_results.assert_called_once_with(
        {"rdf": mock_result.rdf, "rdf_grid": mock_result.rdf_grid},
        conflict_mode=mocker.ANY,
    )
    assert db_handler.return_value.scheme_tables.run_id == 1
    assert read_run.call_args.kwargs["database_name"] is None
    db_handler.assert_called_once_with(None)


@pytest.mark.unit
def test_compute_rdf_passes_custom_grid(mocker, read_run):
    mock_data = mocker.Mock()
    mock_data.inputs = {"dimension": {"value": "_3D"}}
    mock_data.results = {}
    read_run.return_value = mock_data
    mock_result = mocker.Mock()
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_result)
    mocker.patch("qupled.util.dimension.Dimension.from_dict")
    custom_grid = np.array([0.5, 1.0, 1.5])
    compute_rdf(run_id=2, rdf_grid=custom_grid)
    args = mock_result.compute_rdf.call_args[0]
    assert np.allclose(args[1], custom_grid)


@pytest.mark.unit
def test_compute_rdf_passes_database_name(mocker, read_run, db_handler):
    mock_data = mocker.Mock()
    mock_data.inputs = {"dimension": {"value": "_3D"}}
    mock_data.results = {}
    read_run.return_value = mock_data
    mock_result = mocker.Mock()
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_result)
    mocker.patch("qupled.util.dimension.Dimension.from_dict")
    database_name = "test_db"
    compute_rdf(run_id=5, database_name=database_name)
    assert read_run.call_args.kwargs["database_name"] == database_name
    db_handler.assert_called_once_with(database_name)


@pytest.mark.unit
def test_compute_itcf_hf_theory(mocker, read_run, db_handler, scheme_tables):
    mock_data = mocker.Mock()
    mock_data.inputs = {}
    mock_data.results = {}
    mock_data.run = {"theory": "HF"}
    read_run.return_value = mock_data
    mock_results = mocker.Mock()
    mock_results.itcf = np.array([[1.0, 2.0]])
    mock_results.tau = np.array([0.0, 0.1])
    mocker.patch("qupled.schemes.hf.Input.from_dict", return_value=mocker.Mock())
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_results)
    compute_itcf(run_id=1)
    mock_results.compute_itcf.assert_called_once()
    scheme_tables.insert_results.assert_called_once_with(
        {"itcf": mock_results.itcf, "tau": mock_results.tau},
        conflict_mode=mocker.ANY,
    )
    assert read_run.call_args.kwargs["database_name"] is None
    db_handler.assert_called_once_with(None)


@pytest.mark.unit
def test_compute_itcf_rpa_theory(mocker, read_run, scheme_tables):
    mock_data = mocker.Mock()
    mock_data.inputs = {}
    mock_data.results = {}
    mock_data.run = {"theory": "RPA"}
    read_run.return_value = mock_data
    mock_results = mocker.Mock()
    mock_results.itcf = np.array([[3.0, 4.0]])
    mock_results.tau = np.array([0.0, 0.2])
    mocker.patch("qupled.schemes.hf.Input.from_dict", return_value=mocker.Mock())
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_results)
    compute_itcf(run_id=2)
    mock_results.compute_itcf.assert_called_once()
    scheme_tables.insert_results.assert_called_once_with(
        {"itcf": mock_results.itcf, "tau": mock_results.tau},
        conflict_mode=mocker.ANY,
    )


@pytest.mark.unit
def test_compute_itcf_passes_custom_tau(mocker, read_run):
    mock_data = mocker.Mock()
    mock_data.inputs = {}
    mock_data.results = {}
    mock_data.run = {"theory": "HF"}
    read_run.return_value = mock_data
    mock_results = mocker.Mock()
    mocker.patch("qupled.schemes.hf.Input.from_dict", return_value=mocker.Mock())
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_results)
    custom_tau = np.array([0.0, 0.1, 0.2])
    compute_itcf(run_id=3, tau=custom_tau)
    call_args = mock_results.compute_itcf.call_args[0]
    assert np.allclose(call_args[1], custom_tau)


@pytest.mark.unit
def test_compute_itcf_passes_database_name(mocker, read_run, db_handler):
    mock_data = mocker.Mock()
    mock_data.inputs = {}
    mock_data.results = {}
    mock_data.run = {"theory": "HF"}
    read_run.return_value = mock_data
    mock_results = mocker.Mock()
    mocker.patch("qupled.schemes.hf.Input.from_dict", return_value=mocker.Mock())
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_results)
    database_name = "test_db"
    compute_itcf(run_id=6, database_name=database_name)
    assert read_run.call_args.kwargs["database_name"] == database_name
    db_handler.assert_called_once_with(database_name)


@pytest.mark.unit
def test_compute_itcf_stores_results_from_result_object(
    mocker, read_run, scheme_tables
):
    """Verify itcf and tau stored in DB come from results.itcf/tau, not the raw return value."""
    mock_data = mocker.Mock()
    mock_data.inputs = {}
    mock_data.results = {}
    mock_data.run = {"theory": "HF"}
    read_run.return_value = mock_data
    mock_results = mocker.Mock()
    mock_results.compute_itcf.return_value = None  # method returns None
    mock_results.itcf = np.array([[5.0, 6.0]])
    mock_results.tau = np.array([0.0, 0.3])
    mocker.patch("qupled.schemes.hf.Input.from_dict", return_value=mocker.Mock())
    mocker.patch("qupled.schemes.hf.Result.from_dict", return_value=mock_results)
    compute_itcf(run_id=4)
    stored = scheme_tables.insert_results.call_args[0][0]
    assert np.allclose(stored["itcf"], mock_results.itcf)
    assert np.allclose(stored["tau"], mock_results.tau)
