import os
import zipfile
import shutil
import glob
import pytest
import numpy as np
from qupled.native import QVSStls as NativeQVSStls
from qupled.util import Hdf, MPI
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 0


def test_cleanFixedAdrFiles_with_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    mockIsDir = mocker.patch.object(os.path, os.path.isdir.__name__, return_value=True)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockIsDir.call_count == 1
    assert mockRemove.call_count == 1


def test_save(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    try:
        scheme = NativeQVSStls(qvsstls_input.toNative())
        qvsstls.hdfFileName = qvsstls._getHdfFile(scheme.inputs)
        qvsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qvsstls.hdfFileName)
        inspectData = Hdf().inspect(qvsstls.hdfFileName)
        expectedEntries = [
            Hdf.EntryKeys.COUPLING,
            Hdf.EntryKeys.DEGENERACY,
            Hdf.EntryKeys.THEORY,
            Hdf.EntryKeys.ERROR,
            Hdf.EntryKeys.RESOLUTION,
            Hdf.EntryKeys.CUTOFF,
            Hdf.EntryKeys.FREQUENCY_CUTOFF,
            Hdf.EntryKeys.MATSUBARA,
            Hdf.EntryKeys.IDR,
            Hdf.EntryKeys.SDR,
            Hdf.EntryKeys.SLFC,
            Hdf.EntryKeys.SSF,
            Hdf.EntryKeys.SSF_HF,
            Hdf.EntryKeys.WVG,
            Hdf.EntryKeys.FXC_GRID,
            Hdf.EntryKeys.FXCI,
            Hdf.EntryKeys.ADR,
            Hdf.EntryKeys.ALPHA,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qvsstls.hdfFileName)


def test_setFreeEnergyIntegrand(mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mocker.patch.object(
        Hdf,
        Hdf.read.__name__,
        return_value={
            Hdf.EntryKeys.FXC_GRID: arr1D,
            Hdf.EntryKeys.FXCI: arr2D,
            Hdf.EntryKeys.ALPHA: arr1D,
        },
    )
    fxc = QVSStls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxc.grid, arr1D)
    assert np.array_equal(fxc.alpha, arr1D)
    assert np.array_equal(fxc.integrand, arr2D)
