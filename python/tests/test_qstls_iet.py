import os
from shutil import rmtree
import pytest
import numpy as np
import zipfile as zf
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.quantum import QstlsIet, QuantumIterativeScheme


@pytest.fixture
def qstls_iet_instance():
    return QstlsIet(QstlsIet.Input(1.0, 1.0, "QSTLS-HNC"))


def test_default(qstls_iet_instance):
    assert issubclass(QstlsIet, QuantumIterativeScheme)
    assert issubclass(QstlsIet, qp.Qstls)
    assert qstls_iet_instance.fixedIetSourceFile == ""
    assert qstls_iet_instance.hdfFileName == "rs1.000_theta1.000_QSTLS-HNC.h5"


def test_compute(qstls_iet_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockUnpackFixedAdrFiles = mocker.patch(
        "qupled.quantum.QstlsIet._unpackFixedAdrFiles"
    )
    mockCompute = mocker.patch("qupled.qupled.Qstls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.quantum.QstlsIet._checkStatusAndClean"
    )
    mockSave = mocker.patch("qupled.quantum.QstlsIet._save")
    qstls_iet_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockUnpackFixedAdrFiles.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_unpackFixedAdrFiles(qstls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    inputs = qstls_iet_instance.inputs
    inputs.fixediet = "testFile.zip"
    qstls = QstlsIet(inputs)
    try:
        filenames = ["testFile1.txt", "testFile2.txt", "testFile3.txt"]
        for filename in filenames:
            with open(filename, "w") as file:
                file.write("this is a test file\n")
        with zf.ZipFile(qstls.fixedIetSourceFile, "w") as zipFile:
            for filename in filenames:
                zipFile.write(filename)
                os.remove(filename)
        qstls._unpackFixedAdrFiles()
        assert mockMPIIsRoot.call_count == 1
    finally:
        for filename in filenames:
            if os.path.isfile(filename):
                os.remove(filename)
            if os.path.isfile("testFile.zip"):
                os.remove("testFile.zip")
            if os.path.isdir(qstls.inputs.fixediet):
                rmtree(qstls.inputs.fixediet)


def test_checkStatusAndClean(qstls_iet_instance, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qstls_iet_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 2
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        qstls_iet_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"
    inputs = qstls_iet_instance.inputs
    inputs.fixediet = "testFile.zip"
    qstls = QstlsIet(inputs)
    qstls._checkStatusAndClean(0)
    assert not os.path.isdir(qstls.inputs.fixediet)
    with pytest.raises(SystemExit) as excinfo:
        qstls_iet_instance._checkStatusAndClean(1)
    assert not os.path.isdir(qstls_iet_instance.inputs.fixediet)


def test_save(qstls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    adrFileName = "adr_fixed_theta%5.3f_matsubara%d_%s.zip" % (
        qstls_iet_instance.inputs.degeneracy,
        qstls_iet_instance.inputs.matsubara,
        qstls_iet_instance.inputs.theory,
    )
    try:
        qstls_iet_instance._save()
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qstls_iet_instance.hdfFileName)
        inspectData = Hdf().inspect(qstls_iet_instance.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "adr",
            "idr",
            "sdr",
            "slfc",
            "bf",
            "ssf",
            "ssfHF",
            "wvg",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
        assert os.path.isfile(adrFileName)
    finally:
        os.remove(qstls_iet_instance.hdfFileName)
        os.remove(adrFileName)


def test_getInitialGuess(qstls_iet_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((10, 128))
    n = 128
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"wvg": arr1D, "ssf": arr1D, "adr": arr2D, "matsubara": n},
    )
    guess = QstlsIet.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr1D)
    assert np.array_equal(guess.ssf, arr1D)
    assert np.array_equal(guess.adr, arr2D)
    assert np.array_equal(guess.matsubara, n)
