import os
import pytest
import zipfile
import glob
import shutil
from qupled.native import Qstls as NativeQstls
from qupled.util import HDF, MPI
from qupled.quantum import QstlsIet


@pytest.fixture
def qstls_iet():
    return QstlsIet()


@pytest.fixture
def qstls_iet_input():
    return QstlsIet.Input(1.0, 1.0, "QSTLS-HNC")


def test_default(qstls_iet):
    assert qstls_iet.hdf_file_name is None


def test_compute(qstls_iet, qstls_iet_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockUnpack = mocker.patch.object(QstlsIet, QstlsIet._unpackFixedAdrFiles.__name__)
    mockCompute = mocker.patch.object(QstlsIet, QstlsIet._compute.__name__)
    mockSave = mocker.patch.object(QstlsIet, QstlsIet._save.__name__)
    mockZip = mocker.patch.object(QstlsIet, QstlsIet._zipFixedAdrFiles.__name__)
    mockClean = mocker.patch.object(QstlsIet, QstlsIet._cleanFixedAdrFiles.__name__)
    qstls_iet.compute(qstls_iet_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockUnpack.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1
    assert mockZip.call_count == 1
    assert mockClean.call_count == 1


def test_unpackFixedAdrFiles_no_files(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockExtractAll = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qstls_iet._unpackFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 0
    assert mockExtractAll.call_count == 0


def test_unpackFixedAdrFiles_with_files(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockExtractAll = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qstls_iet_input.fixed_iet = "testFile.zip"
    qstls_iet._unpackFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockExtractAll.call_count == 1


def test_zipFixedAdrFiles_no_file(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockZip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mockGlob = mocker.patch.object(
        glob, glob.glob.__name__, return_value={"binFile1", "binFile2"}
    )
    mockRemove = mocker.patch.object(os, os.remove.__name__)
    mockWrite = mocker.patch.object(zipfile.ZipFile, zipfile.ZipFile.write.__name__)
    qstls_iet._zipFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockGlob.call_count == 1
    assert mockRemove.call_count == 2
    assert mockWrite.call_count == 2


def test_cleanFixedAdrFiles_no_files(
    qstls_iet,
    qstls_iet_input,
    mocker,
):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qstls_iet._cleanFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 0


def test_cleanFixedAdrFiles_with_files(
    qstls_iet,
    qstls_iet_input,
    mocker,
):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockIsDir = mocker.patch.object(os.path, os.path.isdir.__name__, return_value=True)
    mockRemove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qstls_iet._cleanFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockIsDir.call_count == 1
    assert mockRemove.call_count == 1


def test_save(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQstls(qstls_iet_input.toNative())
        qstls_iet.hdf_file_name = qstls_iet._get_hdf_file(scheme.inputs)
        qstls_iet._save(scheme)
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qstls_iet.hdf_file_name)
        inspectData = HDF().inspect(qstls_iet.hdf_file_name)
        expectedEntries = [
            HDF.EntryKeys.COUPLING.value,
            HDF.EntryKeys.DEGENERACY.value,
            HDF.EntryKeys.THEORY.value,
            HDF.EntryKeys.ERROR.value,
            HDF.EntryKeys.RESOLUTION.value,
            HDF.EntryKeys.CUTOFF.value,
            HDF.EntryKeys.FREQUENCY_CUTOFF.value,
            HDF.EntryKeys.MATSUBARA.value,
            HDF.EntryKeys.ADR.value,
            HDF.EntryKeys.IDR.value,
            HDF.EntryKeys.SDR.value,
            HDF.EntryKeys.SLFC.value,
            HDF.EntryKeys.BF.value,
            HDF.EntryKeys.SSF.value,
            HDF.EntryKeys.SSF_HF.value,
            HDF.EntryKeys.WVG.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qstls_iet.hdf_file_name)
