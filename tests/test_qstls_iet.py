import os
import pytest
import zipfile
import glob
import shutil
from qupled.native import Qstls as NativeQstls
from qupled.util import HDF, MPI
from qupled.qstlsiet import QstlsIet


@pytest.fixture
def qstls_iet():
    return QstlsIet()


@pytest.fixture
def qstls_iet_input():
    return QstlsIet.Input(1.0, 1.0, "QSTLS-HNC")


def test_default(qstls_iet):
    assert qstls_iet.hdf_file_name is None


def test_compute(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_unpack = mocker.patch.object(
        QstlsIet, QstlsIet._unpack_fixed_adr_files.__name__
    )
    mock_compute = mocker.patch.object(QstlsIet, QstlsIet._compute.__name__)
    mock_save = mocker.patch.object(QstlsIet, QstlsIet._save.__name__)
    mock_zip = mocker.patch.object(QstlsIet, QstlsIet._zip_fixed_adr_files.__name__)
    mock_clean = mocker.patch.object(QstlsIet, QstlsIet._clean_fixed_adr_files.__name__)
    qstls_iet.compute(qstls_iet_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_unpack.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_clean.call_count == 1


def test_unpack_fixed_adr_files_no_files(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_extract_all = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qstls_iet._unpack_fixed_adr_files(qstls_iet_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 0
    assert mock_extract_all.call_count == 0


def test_unpack_fixed_adr_files_with_files(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_extract_all = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qstls_iet_input.fixed_iet = "testFile.zip"
    qstls_iet._unpack_fixed_adr_files(qstls_iet_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_extract_all.call_count == 1


def test_zip_fixed_adr_files_no_file(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_glob = mocker.patch.object(
        glob, glob.glob.__name__, return_value={"binFile1", "binFile2"}
    )
    mock_remove = mocker.patch.object(os, os.remove.__name__)
    mock_write = mocker.patch.object(zipfile.ZipFile, zipfile.ZipFile.write.__name__)
    qstls_iet._zip_fixed_adr_files(qstls_iet_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_glob.call_count == 1
    assert mock_remove.call_count == 2
    assert mock_write.call_count == 2


def test_clean_fixed_adr_files_no_files(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_remove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qstls_iet._clean_fixed_adr_files(qstls_iet_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_remove.call_count == 0


def test_clean_fixed_adr_files_with_files(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_is_dir = mocker.patch.object(
        os.path, os.path.isdir.__name__, return_value=True
    )
    mock_remove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qstls_iet._clean_fixed_adr_files(qstls_iet_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_is_dir.call_count == 1
    assert mock_remove.call_count == 1


def test_save(qstls_iet, qstls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQstls(qstls_iet_input.to_native())
        qstls_iet.hdf_file_name = qstls_iet._get_hdf_file(scheme.inputs)
        qstls_iet._save(scheme)
        assert mock_mpi_is_root.call_count == 4
        assert os.path.isfile(qstls_iet.hdf_file_name)
        inspect_data = HDF().inspect(qstls_iet.hdf_file_name)
        expected_entries = [
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
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(qstls_iet.hdf_file_name)
