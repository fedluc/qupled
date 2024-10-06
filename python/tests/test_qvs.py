import os
from shutil import rmtree
import pytest
import numpy as np
import zipfile as zf
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.quantum import QVSStls, QuantumIterativeScheme


@pytest.fixture
def qvsstls_instance():
    return QVSStls(QVSStls.Input(1.0, 1.0))


def test_default(qvsstls_instance):
    assert issubclass(QVSStls, QuantumIterativeScheme)
    assert issubclass(QVSStls, qp.QVSStls)
    assert qvsstls_instance.hdfFileName == "rs1.000_theta1.000_QVSSTLS.h5"

    
def test_compute(qvsstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.QVSStls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.quantum.QVSStls._checkStatusAndClean"
    )
    mockSave = mocker.patch("qupled.quantum.QVSStls._save")
    qvsstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1
    

def test_unpackFixedAdrFiles(qvsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    inputs = qvsstls_instance.inputs
    inputs.fixed = "testFile.zip"
    qvsstls = QVSStls(inputs)
    try:
        filenames = ["testFile1.txt", "testFile2.txt", "testFile3.txt"]
        for filename in filenames:
            with open(filename, "w") as file:
                file.write("this is a test file\n")
        with zf.ZipFile(qvsstls.inputs.fixed, "w") as zipFile:
            for filename in filenames:
                zipFile.write(filename)
                os.remove(filename)
        qvsstls._unpackFixedAdrFiles()
        assert mockMPIIsRoot.call_count == 1
    finally:
        for filename in filenames:
            if os.path.isfile(filename):
                os.remove(filename)
            if os.path.isfile("testFile.zip"):
                os.remove("testFile.zip")
            if os.path.isdir(qvsstls.inputs.fixed):
                rmtree(qvsstls.inputs.tmpRunDir)


def test_checkStatusAndClean(qvsstls_instance, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qvsstls_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 2
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        qvsstls_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"
    inputs = qvsstls_instance.inputs
    inputs.fixed = "testFile.zip"
    qvsstls = QVSStls(inputs)
    qvsstls._checkStatusAndClean(0)
    assert not os.path.isdir(qvsstls_instance.inputs.fixed)
    with pytest.raises(SystemExit) as excinfo:
        qvsstls_instance._checkStatusAndClean(1)
    assert not os.path.isdir(qvsstls_instance.inputs.fixed)


def test_save(qvsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    adrFileName = "adr_fixed_theta%5.3f_matsubara%d.zip" % (
        qvsstls_instance.inputs.degeneracy,
        qvsstls_instance.inputs.matsubara,
    )
    try:
        qvsstls_instance._save()
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qvsstls_instance.hdfFileName)
        inspectData = Hdf().inspect(qvsstls_instance.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "idr",
            "sdr",
            "slfc",
            "ssf",
            "ssfHF",
            "wvg",
            "fxcGrid",
            "fxci",
            "adr",
            "alpha",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
        assert os.path.isfile(adrFileName)
    finally:
        os.remove(qvsstls_instance.hdfFileName)
        os.remove(adrFileName)


def test_setFreeEnergyIntegrand(qvsstls_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"fxcGrid": arr1D, "fxci": arr2D, "alpha": arr1D},
    )
    fxc = QVSStls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxc.grid, arr1D)
    assert np.array_equal(fxc.alpha, arr1D)
    assert np.array_equal(fxc.integrand, arr2D)
