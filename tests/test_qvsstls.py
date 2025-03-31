import pytest
import os
import zipfile
import qupled.qstls as qstls
import qupled.qvsstls as qvsstls
import qupled.vsstls as vsstls


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qvsstls.QVSStls()


def test_qvsstls_inheritance():
    assert issubclass(qvsstls.QVSStls, vsstls.VSStls)


def test_compute_with_valid_input(mocker, scheme):
    input = mocker.ANY
    unpack = mocker.patch.object(scheme, "_unpack_fixed_adr_files")
    zip = mocker.patch.object(scheme, "_zip_fixed_adr_files")
    clean = mocker.patch.object(scheme, "_clean_fixed_adr_files")
    super_compute = mocker.patch("qupled.vsstls.VSStls.compute")
    scheme.compute(input)
    unpack.assert_called_once_with(input)
    super_compute.assert_called_once_with(input)
    zip.assert_called_once_with(input)
    clean.assert_called_once_with(input)


def test_compute_handles_exceptions(mocker, scheme, input):
    unpack = mocker.patch.object(scheme, "_unpack_fixed_adr_files")
    zip = mocker.patch.object(scheme, "_zip_fixed_adr_files")
    clean = mocker.patch.object(scheme, "_clean_fixed_adr_files")
    super_compute = mocker.patch(
        "qupled.vsstls.VSStls.compute",
        side_effect=RuntimeError("Test exception"),
    )
    with pytest.raises(RuntimeError, match="Test exception"):
        scheme.compute(input)
    unpack.assert_called_once_with(input)
    super_compute.assert_called_once_with(input)
    zip.assert_not_called()
    clean.assert_not_called()


def test_unpack_fixed_adr_files(scheme, input):
    input.fixed = "test_fixed.zip"
    with zipfile.ZipFile(input.fixed, "w") as zip_file:
        zip_file.writestr("test_file.txt", "test content")
    scheme._unpack_fixed_adr_files(input)
    assert os.path.isdir(input.fixed)
    assert os.path.isfile(os.path.join(input.fixed, "test_file.txt"))


def test_zip_fixed_adr_files_creates_zip_file(mocker, scheme, input):
    input.fixed = ""
    input.degeneracy = 1.23
    input.matsubara = 5
    input.theory = "QSTLS-HNC"
    adr_file_bin = "THETA_test.bin"
    with open(adr_file_bin, "w") as f:
        f.write("test content")
    mock_glob = mocker.patch("glob.glob", return_value=[adr_file_bin])
    scheme._zip_fixed_adr_files(input)
    adr_file_zip = "adr_fixed_theta1.230_matsubara5_QSTLS-HNC.zip"
    assert os.path.isfile(adr_file_zip)
    with zipfile.ZipFile(adr_file_zip, "r") as zip_file:
        assert adr_file_bin in zip_file.namelist()
    assert not os.path.isfile(adr_file_bin)
    os.remove(adr_file_zip)


def test_clean_fixed_adr_files_removes_directory(mocker, scheme, input):
    input.fixed = "test_directory"
    mock_isdir = mocker.patch("os.path.isdir", return_value=True)
    mock_rmtree = mocker.patch("shutil.rmtree")
    scheme._clean_fixed_adr_files(input)
    mock_isdir.assert_called_once_with(input.fixed)
    mock_rmtree.assert_called_once_with(input.fixed)


def test_clean_fixed_adr_files_no_directory(mocker, scheme, input):
    input.fixed = "non_existent_directory"
    mock_isdir = mocker.patch("os.path.isdir", return_value=False)
    mock_rmtree = mocker.patch("shutil.rmtree")
    scheme._clean_fixed_adr_files(input)
    mock_isdir.assert_called_once_with(input.fixed)
    mock_rmtree.assert_not_called()


def test_get_free_energy_integrand(mocker):
    run_id = mocker.ANY
    database_name = mocker.ANY
    get_free_energy_integrand = mocker.patch(
        "qupled.vsstls.VSStls.get_free_energy_integrand"
    )
    result = qvsstls.QVSStls.get_free_energy_integrand(run_id, database_name)
    get_free_energy_integrand.assert_called_once_with(run_id, database_name)
    assert result == get_free_energy_integrand.return_value


def test_qvsstls_input_inheritance():
    assert issubclass(qvsstls.Input, (qstls.Input, vsstls.Input))


def test_qvsstls_input_initialization_valid_theory(mocker):
    qstls_init = mocker.patch("qupled.qstls.Input.__init__")
    vsstls_init = mocker.patch("qupled.vsstls.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    input = qvsstls.Input(coupling, degeneracy)
    qstls_init.assert_called_once_with(input, coupling, degeneracy)
    vsstls_init.assert_called_once_with(input, coupling, degeneracy)
    assert input.theory == "QVSSTLS"


def test_qvsstls_result_inheritance():
    assert issubclass(qvsstls.Result, (qstls.Result, vsstls.Result))


def test_qvsstls_result_initialization():
    result = qvsstls.Result()
    assert isinstance(result, qvsstls.Result)
