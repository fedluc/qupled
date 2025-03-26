import os
import zipfile
import shutil
import glob
import pytest
import numpy as np
from qupled.native import QVSStls as NativeQVSStls
from qupled.util import DataBase, MPI
from qupled.qvsstls import QVSStls


@pytest.fixture
def qvsstls():
    return QVSStls()


@pytest.fixture
def qvsstls_input():
    return QVSStls.Input(1.0, 1.0)


def test_default(qvsstls):
    assert qvsstls.hdf_file_name is None


def test_compute(qvsstls, qvsstls_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_unpack = mocker.patch.object(QVSStls, QVSStls._unpack_fixed_adr_files.__name__)
    mock_compute = mocker.patch.object(QVSStls, QVSStls._compute.__name__)
    mock_save = mocker.patch.object(QVSStls, QVSStls._save.__name__)
    mock_zip = mocker.patch.object(QVSStls, QVSStls._zip_fixed_adr_files.__name__)
    mock_clean = mocker.patch.object(QVSStls, QVSStls._clean_fixed_adr_files.__name__)
    qvsstls.compute(qvsstls_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_unpack.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_clean.call_count == 1


def test_unpack_fixed_adr_files_no_files(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_extract_all = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qvsstls._unpack_fixed_adr_files(qvsstls_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 0
    assert mock_extract_all.call_count == 0


def test_unpack_fixed_adr_files_with_files(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_extract_all = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.extractall.__name__, return_value=None
    )
    qvsstls_input.fixed = "testFile.zip"
    qvsstls._unpack_fixed_adr_files(qvsstls_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_extract_all.call_count == 1


def test_zip_fixed_adr_files_no_file(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_zip = mocker.patch.object(
        zipfile.ZipFile, zipfile.ZipFile.__init__.__name__, return_value=None
    )
    mock_glob = mocker.patch.object(
        glob, glob.glob.__name__, return_value={"binFile1", "binFile2"}
    )
    mock_remove = mocker.patch.object(os, os.remove.__name__)
    mock_write = mocker.patch.object(zipfile.ZipFile, zipfile.ZipFile.write.__name__)
    qvsstls._zip_fixed_adr_files(qvsstls_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_zip.call_count == 1
    assert mock_glob.call_count == 1
    assert mock_remove.call_count == 2
    assert mock_write.call_count == 2


def test_clean_fixed_adr_files_no_files(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_remove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._clean_fixed_adr_files(qvsstls_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_remove.call_count == 0


def test_clean_fixed_adr_files_with_files(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_is_dir = mocker.patch.object(
        os.path, os.path.isdir.__name__, return_value=True
    )
    mock_remove = mocker.patch.object(shutil, shutil.rmtree.__name__)
    qvsstls._clean_fixed_adr_files(qvsstls_input)
    assert mock_mpi_is_root.call_count == 1
    assert mock_is_dir.call_count == 1
    assert mock_remove.call_count == 1


def test_save(qvsstls, qvsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQVSStls(qvsstls_input.to_native())
        qvsstls.hdf_file_name = qvsstls._get_hdf_file(scheme.inputs)
        qvsstls._save(scheme)
        assert mock_mpi_is_root.call_count == 4
        assert os.path.isfile(qvsstls.hdf_file_name)
        inspect_data = DataBase().inspect(qvsstls.hdf_file_name)
        expected_entries = [
            DataBase.ResultNames.COUPLING.value,
            DataBase.ResultNames.DEGENERACY.value,
            DataBase.ResultNames.THEORY.value,
            DataBase.ResultNames.ERROR.value,
            DataBase.ResultNames.RESOLUTION.value,
            DataBase.ResultNames.CUTOFF.value,
            DataBase.ResultNames.FREQUENCY_CUTOFF.value,
            DataBase.ResultNames.MATSUBARA.value,
            DataBase.ResultNames.IDR.value,
            DataBase.ResultNames.SDR.value,
            DataBase.ResultNames.SLFC.value,
            DataBase.ResultNames.SSF.value,
            DataBase.ResultNames.SSF_HF.value,
            DataBase.ResultNames.WVG.value,
            DataBase.ResultNames.FXC_GRID.value,
            DataBase.ResultNames.FXC_INT.value,
            DataBase.ResultNames.ADR.value,
            DataBase.ResultNames.ALPHA.value,
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(qvsstls.hdf_file_name)


def test_set_free_energy_integrand(mocker):
    arr1d = np.ones(10)
    arr2d = np.ones((3, 10))
    mocker.patch.object(
        DataBase,
        DataBase.read.__name__,
        return_value={
            DataBase.ResultNames.FXC_GRID.value: arr1d,
            DataBase.ResultNames.FXC_INT.value: arr2d,
            DataBase.ResultNames.ALPHA.value: arr1d,
        },
    )
    fxc = QVSStls.get_free_energy_integrand("dummyFileName")
    assert np.array_equal(fxc.grid, arr1d)
    assert np.array_equal(fxc.alpha, arr1d)
    assert np.array_equal(fxc.integrand, arr2d)
