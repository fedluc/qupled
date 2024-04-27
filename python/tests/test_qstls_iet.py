import os
from shutil import rmtree
import pytest
import numpy as np
import zipfile as zf
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.quantum import QstlsIet
from qupled.quantum import Qstls


@pytest.fixture
def qstls_iet_instance():
    return QstlsIet(1.0, 1.0, "QSTLS-HNC")


def test_default(qstls_iet_instance):
    assert issubclass(QstlsIet, Qstls)
    assert all(
        x == y
        for x, y in zip(
            qstls_iet_instance.allowedTheories, ["QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"]
        )
    )
    assert qstls_iet_instance.inputs.coupling == 1.0
    assert qstls_iet_instance.inputs.degeneracy == 1.0
    assert qstls_iet_instance.inputs.theory == "QSTLS-HNC"
    assert all(
        x == y
        for x, y in zip(qstls_iet_instance.inputs.chemicalPotential, [-10.0, 10.0])
    )
    assert qstls_iet_instance.inputs.cutoff == 10.0
    assert qstls_iet_instance.inputs.error == 1.0e-5
    assert qstls_iet_instance.inputs.fixed == ""
    assert qstls_iet_instance.inputs.fixediet == ""
    assert qstls_iet_instance.inputs.iet == "standard"
    assert qstls_iet_instance.inputs.mixing == 1.0
    assert qstls_iet_instance.inputs.iterations == 1000
    assert qstls_iet_instance.inputs.matsubara == 128
    assert qstls_iet_instance.inputs.outputFrequency == 10
    assert qstls_iet_instance.inputs.recoveryFile == ""
    assert qstls_iet_instance.inputs.resolution == 0.1
    assert qstls_iet_instance.inputs.int2DScheme == "full"
    assert qstls_iet_instance.inputs.threads == 1
    assert qstls_iet_instance.inputs.intError == 1.0e-5
    assert qstls_iet_instance.tmpRunDir is None
    assert qstls_iet_instance.scheme is None
    assert qstls_iet_instance.hdfFileName is None


def test_set_input():
    qstls = QstlsIet(
        2.0,
        0.5,
        "QSTLS-LCT",
        [-5, 5],
        20,
        1.0e-8,
        "fixedAdrFile",
        "fixedAdrIetFile",
        "sqrt",
        0.5,
        None,
        100,
        32,
        100,
        "recoveryFile",
        0.01,
        "segregated",
        16,
    )
    assert qstls.inputs.coupling == 2.0
    assert qstls.inputs.degeneracy == 0.5
    assert qstls.inputs.theory == "QSTLS-LCT"
    assert all(x == y for x, y in zip(qstls.inputs.chemicalPotential, [-5.0, 5.0]))
    assert qstls.inputs.cutoff == 20
    assert qstls.inputs.error == 1.0e-8
    assert qstls.inputs.fixed == "fixedAdrFile"
    assert qstls.inputs.fixediet == "fixedAdrIetFile"
    assert qstls.inputs.iet == "sqrt"
    assert qstls.inputs.mixing == 0.5
    assert qstls.inputs.iterations == 100
    assert qstls.inputs.matsubara == 32
    assert qstls.inputs.outputFrequency == 100
    assert qstls.inputs.recoveryFile == "recoveryFile"
    assert qstls.inputs.resolution == 0.01
    assert qstls.inputs.int2DScheme == "segregated"
    assert qstls.inputs.threads == 16


def test_compute(qstls_iet_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.quantum.QstlsIet._checkInputs")
    mockSetFixedIetFileName = mocker.patch(
        "qupled.quantum.QstlsIet._setFixedIetFileName"
    )
    mockCompute = mocker.patch("qupled.qupled.Qstls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.quantum.QstlsIet._checkStatusAndClean"
    )
    mockSetHdfFile = mocker.patch("qupled.quantum.QstlsIet._setHdfFile")
    mockSave = mocker.patch("qupled.quantum.QstlsIet._save")
    qstls_iet_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockSetFixedIetFileName.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_setFixedIetFileName(qstls_iet_instance, mocker):
    mockUnpackFixedAdrFiles = mocker.patch(
        "qupled.quantum.QstlsIet._unpackFixedAdrFiles"
    )
    qstls_iet_instance._setFixedIetFileName()
    assert qstls_iet_instance.tmpRunDir is None
    assert qstls_iet_instance.inputs.fixediet == ""
    assert mockUnpackFixedAdrFiles.call_count == 0
    qstls_iet_instance.inputs.fixediet = "fixedAdrIetFile"
    qstls_iet_instance._setFixedIetFileName()
    assert qstls_iet_instance.tmpRunDir == "qupled_tmp_run_directory"
    assert qstls_iet_instance.inputs.fixediet == qstls_iet_instance.tmpRunDir
    assert mockUnpackFixedAdrFiles.call_count == 1


def test_unpackFixedAdrFiles(qstls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qstls_iet_instance.inputs.fixediet = "testFile.zip"
    qstls_iet_instance.tmpRunDir = "qupled_tmp_run_directory"
    try:
        filenames = ["testFile1.txt", "testFile2.txt", "testFile3.txt"]
        for filename in filenames:
            with open(filename, "w") as file:
                file.write("this is a test file\n")
        with zf.ZipFile(qstls_iet_instance.inputs.fixediet, "w") as zipFile:
            for filename in filenames:
                zipFile.write(filename)
                os.remove(filename)
        qstls_iet_instance._unpackFixedAdrFiles()
        assert mockMPIIsRoot.call_count == 1
    finally:
        for filename in filenames:
            if os.path.isfile(filename):
                os.remove(filename)
            if os.path.isfile("testFile.zip"):
                os.remove("testFile.zip")
            if os.path.isdir(qstls_iet_instance.tmpRunDir):
                rmtree(qstls_iet_instance.tmpRunDir)


def test_checkStatusAndClean(qstls_iet_instance, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qstls_iet_instance.scheme = qp.Qstls(qstls_iet_instance.inputs)
    qstls_iet_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 2
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        qstls_iet_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"
    qstls_iet_instance.tmpRunDir = "tmpRunDir"
    os.makedirs(qstls_iet_instance.tmpRunDir, exist_ok=True)
    assert os.path.isdir(qstls_iet_instance.tmpRunDir)
    qstls_iet_instance._checkStatusAndClean(0)
    assert not os.path.isdir(qstls_iet_instance.tmpRunDir)
    os.makedirs(qstls_iet_instance.tmpRunDir, exist_ok=True)
    assert os.path.isdir(qstls_iet_instance.tmpRunDir)
    with pytest.raises(SystemExit) as excinfo:
        qstls_iet_instance._checkStatusAndClean(1)
    assert not os.path.isdir(qstls_iet_instance.tmpRunDir)


def test_save(qstls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qstls_iet_instance.scheme = qp.Qstls(qstls_iet_instance.inputs)
    qstls_iet_instance._setHdfFile()
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


def test_setGuess(qstls_iet_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((10, 128))
    n = 128
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"wvg": arr1D, "ssf": arr1D, "adr": arr2D, "matsubara": n},
    )
    qstls_iet_instance.setGuess("dummyFileName")
    assert np.array_equal(qstls_iet_instance.inputs.guess.wvg, arr1D)
    assert np.array_equal(qstls_iet_instance.inputs.guess.ssf, arr1D)
    assert np.array_equal(qstls_iet_instance.inputs.guess.adr, arr2D)
    assert np.array_equal(qstls_iet_instance.inputs.guess.matsubara, n)
