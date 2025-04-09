import os
import zipfile

import pytest

import qupled.qstls as qstls
import qupled.qstlsiet as qstlsiet
import qupled.stlsiet as stlsiet


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qstlsiet.QstlsIet()


def test_qstls_iet_inheritance():
    assert issubclass(qstlsiet.QstlsIet, qstls.Qstls)


def test_qstls_iet_initialization(mocker):
    super_init = mocker.patch("qupled.qstls.Qstls.__init__")
    scheme = qstlsiet.QstlsIet()
    super_init.assert_called_once()
    assert isinstance(scheme.results, qstlsiet.Result)


def test_compute_with_valid_input(mocker, scheme):
    input = mocker.ANY
    unpack = mocker.patch.object(scheme, "_unpack_fixed_adr_files")
    zip = mocker.patch.object(scheme, "_zip_fixed_adr_files")
    clean = mocker.patch.object(scheme, "_clean_fixed_adr_files")
    super_compute = mocker.patch("qupled.qstls.Qstls.compute")
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
        "qupled.qstls.Qstls.compute",
        side_effect=RuntimeError("Test exception"),
    )
    with pytest.raises(RuntimeError, match="Test exception"):
        scheme.compute(input)
    unpack.assert_called_once_with(input)
    super_compute.assert_called_once_with(input)
    zip.assert_not_called()
    clean.assert_not_called()


def test_unpack_fixed_adr_files(scheme, input):
    input.fixed_iet = "test_fixed_iet.zip"
    with zipfile.ZipFile(input.fixed_iet, "w") as zip_file:
        zip_file.writestr("test_file.txt", "test content")
    scheme._unpack_fixed_adr_files(input)
    assert os.path.isdir(input.fixed_iet)
    assert os.path.isfile(os.path.join(input.fixed_iet, "test_file.txt"))


def test_zip_fixed_adr_files(scheme, input):
    input.fixed_iet = ""
    input.degeneracy = 1.23
    input.matsubara = 5
    input.theory = "QSTLS-HNC"
    adr_file_bin = "adr_fixed_theta1.230_matsubara5_QSTLS-HNC_wv1.bin"
    with open(adr_file_bin, "w") as f:
        f.write("test content")
    scheme._zip_fixed_adr_files(input)
    adr_file_zip = "adr_fixed_theta1.230_matsubara5_QSTLS-HNC.zip"
    assert os.path.isfile(adr_file_zip)
    with zipfile.ZipFile(adr_file_zip, "r") as zip_file:
        assert (
            "adr_fixed_theta1.230_matsubara5_QSTLS-HNC_wv1.bin" in zip_file.namelist()
        )


def test_clean_fixed_adr_files_removes_directory(mocker, scheme, input):
    input.fixed_iet = "test_directory"
    mock_isdir = mocker.patch("os.path.isdir", return_value=True)
    mock_rmtree = mocker.patch("shutil.rmtree")
    scheme._clean_fixed_adr_files(input)
    mock_isdir.assert_called_once_with(input.fixed_iet)
    mock_rmtree.assert_called_once_with(input.fixed_iet)


def test_clean_fixed_adr_files_no_directory(mocker, scheme, input):
    input.fixed_iet = "non_existent_directory"
    mock_isdir = mocker.patch("os.path.isdir", return_value=False)
    mock_rmtree = mocker.patch("shutil.rmtree")
    scheme._clean_fixed_adr_files(input)
    mock_isdir.assert_called_once_with(input.fixed_iet)
    mock_rmtree.assert_not_called()


def test_qstls_iet_input_inheritance():
    assert issubclass(qstlsiet.Input, (stlsiet.Input, qstls.Input))


def test_qstls_iet_input_initialization_valid_theory(mocker):
    qstls_init = mocker.patch("qupled.qstls.Input.__init__")
    stls_iet_init = mocker.patch("qupled.stlsiet.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    theory = "QSTLS-HNC"
    input = qstlsiet.Input(coupling, degeneracy, theory)
    qstls_init.assert_called_once_with(input, coupling, degeneracy)
    stls_iet_init.assert_called_once_with(input, coupling, degeneracy, "STLS-HNC")
    assert input.theory == theory
    assert input.fixed_iet == ""


def test_qstls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        qstlsiet.Input(1.0, 1.0, "INVALID-THEORY")


def test_qstls_iet_result_inheritance():
    assert issubclass(qstlsiet.Result, (stlsiet.Result, qstls.Result))


def test_qstls_iet_result_initialization():
    result = qstlsiet.Result()
    assert isinstance(result, qstlsiet.Result)
