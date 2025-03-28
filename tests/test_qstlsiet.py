import pytest
import os
import zipfile
from qupled.base import QuantumIterativeScheme
from qupled.native import Qstls as NativeQstls
from qupled.qstls import Qstls
from qupled.qstlsiet import QstlsIet
from qupled.stlsiet import StlsIet


@pytest.fixture
def qstls_iet_input():
    return QstlsIet.Input(coupling=1.0, degeneracy=2.0, theory="QSTLS-HNC")


@pytest.fixture
def qstls_iet_instance():
    return QstlsIet()


def test_qstls_iet_inheritance():
    assert issubclass(QstlsIet, QuantumIterativeScheme)


def test_compute_with_valid_input(mocker, qstls_iet_instance, qstls_iet_input):
    native_input = mocker.MagicMock()
    result = mocker.MagicMock()
    unpack = mocker.patch.object(qstls_iet_instance, "_unpack_fixed_adr_files")
    zip = mocker.patch.object(qstls_iet_instance, "_zip_fixed_adr_files")
    clean = mocker.patch.object(qstls_iet_instance, "_clean_fixed_adr_files")
    mocker.patch.object(qstls_iet_instance, "Result", return_value=result)
    QstlsInput = mocker.patch("qupled.native.QstlsInput", return_value=native_input)
    super_compute = mocker.patch("qupled.base.QuantumIterativeScheme.compute")
    qstls_iet_instance.compute(qstls_iet_input)
    unpack.assert_called_once_with(qstls_iet_input)
    super_compute.assert_called_once_with(
        qstls_iet_input, NativeQstls, native_input, result
    )
    zip.assert_called_once_with(qstls_iet_input)
    clean.assert_called_once_with(qstls_iet_input)


def test_compute_handles_exceptions(mocker, qstls_iet_instance, qstls_iet_input):
    native_input = mocker.MagicMock()
    result = mocker.MagicMock()
    unpack = mocker.patch.object(qstls_iet_instance, "_unpack_fixed_adr_files")
    zip = mocker.patch.object(qstls_iet_instance, "_zip_fixed_adr_files")
    clean = mocker.patch.object(qstls_iet_instance, "_clean_fixed_adr_files")
    mocker.patch.object(qstls_iet_instance, "Result", return_value=result)
    mocker.patch("qupled.native.QstlsInput", return_value=native_input)
    super_compute = mocker.patch(
        "qupled.base.QuantumIterativeScheme.compute",
        side_effect=RuntimeError("Test exception"),
    )
    with pytest.raises(RuntimeError, match="Test exception"):
        qstls_iet_instance.compute(qstls_iet_input)
    unpack.assert_called_once_with(qstls_iet_input)
    super_compute.assert_called_once_with(
        qstls_iet_input, mocker.ANY, mocker.ANY, mocker.ANY
    )
    zip.assert_not_called()
    clean.assert_not_called()


def test_unpack_fixed_adr_files(qstls_iet_instance, qstls_iet_input):
    qstls_iet_input.fixed_iet = "test_fixed_iet.zip"
    with zipfile.ZipFile(qstls_iet_input.fixed_iet, "w") as zip_file:
        zip_file.writestr("test_file.txt", "test content")
    qstls_iet_instance._unpack_fixed_adr_files(qstls_iet_input)
    assert os.path.isdir(qstls_iet_input.fixed_iet)
    assert os.path.isfile(os.path.join(qstls_iet_input.fixed_iet, "test_file.txt"))


def test_zip_fixed_adr_files(qstls_iet_instance, qstls_iet_input):
    qstls_iet_input.fixed_iet = ""
    qstls_iet_input.degeneracy = 1.23
    qstls_iet_input.matsubara = 5
    qstls_iet_input.theory = "QSTLS-HNC"
    adr_file_bin = "adr_fixed_theta1.230_matsubara5_QSTLS-HNC_wv1.bin"
    with open(adr_file_bin, "w") as f:
        f.write("test content")
    qstls_iet_instance._zip_fixed_adr_files(qstls_iet_input)
    adr_file_zip = "adr_fixed_theta1.230_matsubara5_QSTLS-HNC.zip"
    assert os.path.isfile(adr_file_zip)
    with zipfile.ZipFile(adr_file_zip, "r") as zip_file:
        assert (
            "adr_fixed_theta1.230_matsubara5_QSTLS-HNC_wv1.bin" in zip_file.namelist()
        )


def test_clean_fixed_adr_files(qstls_iet_instance, qstls_iet_input):
    qstls_iet_input.fixed_iet = "test_directory"
    os.mkdir(qstls_iet_input.fixed_iet)
    qstls_iet_instance._clean_fixed_adr_files(qstls_iet_input)
    assert not os.path.isdir(qstls_iet_input.fixed_iet)


def test_qstls_iet_input_inheritance():
    assert issubclass(QstlsIet.Input, (StlsIet.Input, Qstls.Input))


def test_qstls_iet_input_initialization_valid_theory(mocker):
    qstls_init = mocker.patch("qupled.qstls.Qstls.Input.__init__")
    stls_iet_init = mocker.patch("qupled.stlsiet.StlsIet.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    theory = "QSTLS-HNC"
    valid_input = QstlsIet.Input(coupling, degeneracy, theory)
    qstls_init.assert_called_once_with(valid_input, coupling, degeneracy)
    stls_iet_init.assert_called_once_with(valid_input, coupling, degeneracy, "STLS-HNC")
    assert valid_input.theory == theory
    assert valid_input.fixed_iet == ""


def test_qstls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        QstlsIet.Input(1.0, 1.0, "INVALID-THEORY")


def test_qstls_iet_result_inheritance():
    assert issubclass(QstlsIet.Result, (StlsIet.Result, Qstls.Result))


def test_qstls_iet_result_initialization():
    result = QstlsIet.Result()
    assert isinstance(result, QstlsIet.Result)
