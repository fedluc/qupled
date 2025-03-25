import os

import pytest

from qupled.stlsiet import StlsIet
from qupled.native import Stls as NativeStls
from qupled.util import HDF, MPI


@pytest.fixture
def stls_iet():
    return StlsIet()


@pytest.fixture
def stls_iet_input():
    return StlsIet.Input(1.0, 1.0, "STLS-HNC")


def test_default(stls_iet):
    assert stls_iet.hdf_file_name is None


def test_compute(stls_iet, stls_iet_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(StlsIet, StlsIet._compute.__name__)
    mock_save = mocker.patch.object(StlsIet, StlsIet._save.__name__)
    stls_iet.compute(stls_iet_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1


def test_save(stls_iet, stls_iet_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeStls(stls_iet_input.to_native())
        stls_iet.hdf_file_name = stls_iet._get_hdf_file(scheme.inputs)
        stls_iet._save(scheme)
        assert mock_mpi_is_root.call_count == 3
        assert os.path.isfile(stls_iet.hdf_file_name)
        inspect_data = HDF().inspect(stls_iet.hdf_file_name)
        expected_entries = [
            HDF.ResultNames.COUPLING.value,
            HDF.ResultNames.DEGENERACY.value,
            HDF.ResultNames.THEORY.value,
            HDF.ResultNames.ERROR.value,
            HDF.ResultNames.RESOLUTION.value,
            HDF.ResultNames.CUTOFF.value,
            HDF.ResultNames.FREQUENCY_CUTOFF.value,
            HDF.ResultNames.MATSUBARA.value,
            HDF.ResultNames.BF.value,
            HDF.ResultNames.IDR.value,
            HDF.ResultNames.SDR.value,
            HDF.ResultNames.SLFC.value,
            HDF.ResultNames.SSF.value,
            HDF.ResultNames.SSF_HF.value,
            HDF.ResultNames.WVG.value,
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(stls_iet.hdf_file_name)
