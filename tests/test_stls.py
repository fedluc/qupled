import os

import numpy as np
import pytest

from qupled.stls import Stls
from qupled.native import Stls as NativeStls
from qupled.util import HDF, MPI


@pytest.fixture
def stls():
    return Stls()


@pytest.fixture
def stls_input():
    return Stls.Input(1.0, 1.0)


def test_default(stls):
    assert stls.hdf_file_name is None


def test_compute(stls, stls_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(Stls, Stls._compute.__name__)
    mock_save = mocker.patch.object(Stls, Stls._save.__name__)
    stls.compute(stls_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1


def test_save(stls, stls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeStls(stls_input.to_native())
        stls.hdf_file_name = stls._get_hdf_file(scheme.inputs)
        stls._save(scheme)
        assert mock_mpi_is_root.call_count == 2
        assert os.path.isfile(stls.hdf_file_name)
        inspect_data = HDF().inspect(stls.hdf_file_name)
        expected_entries = [
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
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(stls.hdf_file_name)


def test_get_initial_guess(mocker):
    arr = np.ones(10)
    mocker.patch.object(
        HDF,
        HDF.read.__name__,
        return_value={HDF.EntryKeys.WVG.value: arr, HDF.EntryKeys.SLFC.value: arr},
    )
    guess = Stls.get_initial_guess("dummy_file_name")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.slfc, arr)
