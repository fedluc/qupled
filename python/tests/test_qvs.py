import os
from shutil import rmtree
import pytest
import numpy as np
import zipfile as zf
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import VSStls
from qupled.quantum import Qstls
from qupled.quantum import QVSStls

@pytest.fixture
def qvsstls_instance():
    return QVSStls(1.0, 1.0)


def test_default(qvsstls_instance):
    assert issubclass(QVSStls, VSStls)
    assert issubclass(QVSStls, Qstls)
    assert all(x == y for x, y in zip(qvsstls_instance.allowedTheories, ["QVSSTLS"]))
    assert qvsstls_instance.inputs.coupling == 1.0
    assert qvsstls_instance.inputs.degeneracy == 1.0
    assert qvsstls_instance.inputs.theory == "QVSSTLS"
    assert all(x == y for x, y in zip(qvsstls_instance.inputs.chemicalPotential, [-10.0, 10.0]))
    assert qvsstls_instance.inputs.cutoff == 10.0
    assert qvsstls_instance.inputs.error == 1.0e-5
    assert qvsstls_instance.inputs.fixed == ""
    assert qvsstls_instance.inputs.mixing == 0.5
    assert qvsstls_instance.inputs.iterations == 1000
    assert qvsstls_instance.inputs.matsubara == 128
    assert qvsstls_instance.inputs.outputFrequency == 10
    assert qvsstls_instance.inputs.recoveryFile == ""
    assert qvsstls_instance.inputs.resolution == 0.1
    assert qvsstls_instance.inputs.int2DScheme == "full"
    assert all(x == y for x, y in zip(qvsstls_instance.inputs.alpha, [0.5, 1.0]))
    assert qvsstls_instance.inputs.couplingResolution == 0.1
    assert qvsstls_instance.inputs.degeneracyResolution == 0.1
    assert qvsstls_instance.inputs.errorAlpha == 1.0e-3
    assert qvsstls_instance.inputs.iterationsAlpha == 50
    assert qvsstls_instance.inputs.intError == 1.0e-5
    assert qvsstls_instance.inputs.threads == 1
    assert qvsstls_instance.tmpRunDir is None
    assert qvsstls_instance.scheme is None
    assert qvsstls_instance.hdfFileName is None

                                      
def test_set_input():
    qvsstls = QVSStls(2.0, 0.5, [-5, 5], 20, 1.0e-8, "fixedAdrFile", 0.5, None, 100,
                    32, 100, "recoveryFile", 0.01, "segregated", [0.8, 1.2],
                    0.1, 0.2, 1.0e-5, 100, 1.0e-8, 9)
    assert qvsstls.inputs.coupling == 2.0
    assert qvsstls.inputs.degeneracy == 0.5
    assert qvsstls.inputs.theory == "QVSSTLS"
    assert all(x == y for x, y in zip(qvsstls.inputs.chemicalPotential, [-5.0, 5.0]))
    assert qvsstls.inputs.cutoff == 20
    assert qvsstls.inputs.error == 1.0e-8
    assert qvsstls.inputs.mixing == 0.5
    assert qvsstls.inputs.fixed == "fixedAdrFile"
    assert qvsstls.inputs.iterations == 100
    assert qvsstls.inputs.matsubara == 32
    assert qvsstls.inputs.outputFrequency == 100
    assert qvsstls.inputs.recoveryFile == "recoveryFile"
    assert qvsstls.inputs.resolution == 0.01
    assert qvsstls.inputs.int2DScheme == "segregated"
    assert all(x == y for x, y in zip(qvsstls.inputs.alpha, [0.8, 1.2]))
    assert qvsstls.inputs.couplingResolution == 0.1
    assert qvsstls.inputs.degeneracyResolution == 0.2
    assert qvsstls.inputs.errorAlpha == 1.0e-5
    assert qvsstls.inputs.iterationsAlpha == 100
    assert qvsstls.inputs.intError == 1.0e-8
    assert qvsstls.inputs.threads == 9
    assert qvsstls.scheme is None
    assert qvsstls.hdfFileName is None
    
def test_compute(qvsstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.quantum.QVSStls._checkInputs")
    mockSetFixedAdrFileName = mocker.patch("qupled.quantum.QVSStls._setFixedAdrFileName")
    mockCompute = mocker.patch("qupled.qupled.QVSStls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.quantum.QVSStls._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.quantum.QVSStls._setHdfFile")
    mockSave = mocker.patch("qupled.quantum.QVSStls._save")
    qvsstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockSetFixedAdrFileName.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1

def test_setFixedAdrFileName(qvsstls_instance, mocker):
    mockUnpackFixedAdrFiles = mocker.patch("qupled.quantum.QVSStls._unpackFixedAdrFiles")
    qvsstls_instance._setFixedAdrFileName()
    assert qvsstls_instance.tmpRunDir is None
    assert qvsstls_instance.inputs.fixed == ""
    assert mockUnpackFixedAdrFiles.call_count == 0
    qvsstls_instance.inputs.fixed = "fixedAdrIetFile"
    qvsstls_instance._setFixedAdrFileName()
    assert qvsstls_instance.tmpRunDir == "qupled_tmp_run_directory"
    assert qvsstls_instance.inputs.fixed == qvsstls_instance.tmpRunDir
    assert mockUnpackFixedAdrFiles.call_count == 1

def test_unpackFixedAdrFiles(qvsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qvsstls_instance.inputs.fixed = "testFile.zip"
    qvsstls_instance.tmpRunDir = "qupled_tmp_run_directory"
    try:
        filenames = ["testFile1.txt", "testFile2.txt", "testFile3.txt"]
        for filename in filenames:
            with open(filename, 'w') as file:
                file.write('this is a test file\n')
        with zf.ZipFile(qvsstls_instance.inputs.fixed, "w") as zipFile:
            for filename in filenames:
                zipFile.write(filename)
                os.remove(filename)
        qvsstls_instance._unpackFixedAdrFiles()
        assert mockMPIIsRoot.call_count == 1
    finally:
        for filename in filenames:
            if os.path.isfile(filename) : os.remove(filename)
            if os.path.isfile("testFile.zip") : os.remove("testFile.zip")
            if os.path.isdir(qvsstls_instance.tmpRunDir) : rmtree(qvsstls_instance.tmpRunDir)

def test_checkStatusAndClean(qvsstls_instance, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qvsstls_instance.scheme = qp.Qstls(qvsstls_instance.inputs)
    qvsstls_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 2
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        qvsstls_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"
    qvsstls_instance.tmpRunDir = "tmpRunDir"
    os.makedirs(qvsstls_instance.tmpRunDir, exist_ok=True)
    assert os.path.isdir(qvsstls_instance.tmpRunDir)
    qvsstls_instance._checkStatusAndClean(0)
    assert not os.path.isdir(qvsstls_instance.tmpRunDir)
    os.makedirs(qvsstls_instance.tmpRunDir, exist_ok=True)
    assert os.path.isdir(qvsstls_instance.tmpRunDir)
    with pytest.raises(SystemExit) as excinfo:
        qvsstls_instance._checkStatusAndClean(1)
    assert not os.path.isdir(qvsstls_instance.tmpRunDir)
    
def test_save(qvsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qvsstls_instance.scheme = qp.QVSStls(qvsstls_instance.inputs)
    qvsstls_instance._setHdfFile()
    adrFileName = "adr_fixed_theta%5.3f_matsubara%d.zip" % (qvsstls_instance.inputs.degeneracy,
                                                            qvsstls_instance.inputs.matsubara)
    try:
        qvsstls_instance._save()
        assert mockMPIIsRoot.call_count == 5
        assert os.path.isfile(qvsstls_instance.hdfFileName)
        inspectData = Hdf().inspect(qvsstls_instance.hdfFileName)
        expectedEntries = ["coupling", "degeneracy", "theory", "error",
                           "resolution", "cutoff", "matsubara",
                           "idr", "sdr", "slfc",  "ssf", "ssfHF",
                           "wvg", "fxcGrid", "fxci", "adr", "alpha"]
        for entry in expectedEntries:
            assert entry in inspectData
        assert os.path.isfile(adrFileName)
    finally:
        os.remove(qvsstls_instance.hdfFileName)
        os.remove(adrFileName)


def test_setFreeEnergyIntegrand(qvsstls_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch("qupled.util.Hdf.read", return_value={"fxcGrid" : arr1D, "fxci" : arr2D, "alpha" : arr1D})
    qvsstls_instance.setFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(qvsstls_instance.inputs.freeEnergyIntegrand.grid, arr1D)
    assert np.array_equal(qvsstls_instance.inputs.freeEnergyIntegrand.alpha, arr1D)
    assert np.array_equal(qvsstls_instance.inputs.freeEnergyIntegrand.integrand, arr2D)
