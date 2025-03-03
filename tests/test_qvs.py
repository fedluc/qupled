import os
import zipfile
import shutil
import glob
import pytest
import numpy as np
from qupled.native import QVSStls as NativeQVSStls
from qupled.util import HDF, MPI
from qupled.quantum import QVSStls


@pytest.fixture
def qvsstls():
    return QVSStls()


@pytest.fixture
def qvsstls_input():
    return QVSStls.Input(1.0, 1.0)


def test_default(qvsstls):
    assert qvsstls.hdfFileName is None


def test_compute(qvsstls, qvsstls_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockUnpack = mocker.patch.object(QVSStls, QVSStls._unpackFixedAdrFiles.__name__)
    mockCompute = mocker.patch.object(QVSStls, QVSStls._compute.__name__)
    mockSave = mocker.patch.object(QVSStls, QVSStls._save.__name__)
    mockZip = mocker.patch.object(QVSStls, QVSStls._zipFixedAdrFiles.__name__)
    mockClean = mocker.patch.object(QVSStls, QVSStls._cleanFixedAdrFiles.__name__)
    qvsstls.compute(qvsstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockUnpack.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1
    assert mockZip.call_count == 1
    assert mockClean.call_count == 1


def test_unpackFixedAdrFiles_no_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockExtractAll = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qvsstls._unpackFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 0
    assert mockExtractAll.call_count == 0


def test_unpackFixedAdrFiles_with_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockExtractAll = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qvsstls_input.fixed = "testFile.zip"
    qvsstls._unpackFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockExtractAll.call_count == 1


def test_zipFixedAdrFiles_no_file(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockGlob = mocker.patch.object(
        glob, glob.glob.__name__, return_value={"binFile1", "binFile2"}
    )
    mockRemove = mocker.patch.object(os, os.remove.__name__)
    mockWrite = mocker.patch.object(zipfile.ZipFile, zipfile.ZipFile.write.__name__)
    qvsstls._zipFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockGlob.call_count == 1
    assert mockRemove.call_count == 2
    assert mockWrite.call_count == 2


def test_cleanFixedAdrFiles_no_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 0


def test_cleanFixedAdrFiles_with_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockIsDir = mocker.patch.object(os.path, os.path.isdir.__name__, return_value=True)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockIsDir.call_count == 1
    assert mockRemove.call_count == 1


def test_save(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQVSStls(qvsstls_input.toNative())
        qvsstls.hdfFileName = qvsstls._getHdfFile(scheme.inputs)
        qvsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qvsstls.hdfFileName)
        inspectData = HDF().inspect(qvsstls.hdfFileName)
        expectedEntries = [
            HDF.EntryKeys.COUPLING.value,
            HDF.EntryKeys.DEGENERACY.value,
            HDF.EntryKeys.THEORY.value,
            HDF.EntryKeys.ERROR.value,
            HDF.EntryKeys.RESOLUTION.value,
            HDF.EntryKeys.CUTOFF.value,
            HDF.EntryKeys.FREQUENCY_CUTOFF.value,
            HDF.EntryKeys.MATSUBARA.value,
            HDF.EntryKeys.IDR.value,
            HDF.EntryKeys.SDR.value,
            HDF.EntryKeys.SLFC.value,
            HDF.EntryKeys.SSF.value,
            HDF.EntryKeys.SSF_HF.value,
            HDF.EntryKeys.WVG.value,
            HDF.EntryKeys.FXC_GRID.value,
            HDF.EntryKeys.FXCI.value,
            HDF.EntryKeys.ADR.value,
            HDF.EntryKeys.ALPHA.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qvsstls.hdfFileName)


def test_setFreeEnergyIntegrand(mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mocker.patch.object(
        HDF,
        HDF.read.__name__,
        return_value={
            HDF.EntryKeys.FXC_GRID.value: arr1D,
            HDF.EntryKeys.FXCI.value: arr2D,
            HDF.EntryKeys.ALPHA.value: arr1D,
        },
    )
    fxc = QVSStls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxc.grid, arr1D)
    assert np.array_equal(fxc.alpha, arr1D)
    assert np.array_equal(fxc.integrand, arr2D)
