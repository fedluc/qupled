import os

import numpy as np
import pytest

from qupled import native
from qupled.vsstls import VSStls
from qupled.util import DataBase, MPI


@pytest.fixture
def vsstls():
    return VSStls()


@pytest.fixture
def vsstls_input():
    return VSStls.Input(1.0, 1.0)


def test_default(vsstls):
    assert vsstls.hdf_file_name is None


def test_compute(vsstls, vsstls_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(VSStls, VSStls._compute.__name__)
    mock_save = mocker.patch.object(VSStls, VSStls._save.__name__)
    vsstls.compute(vsstls_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1


def test_save(vsstls, vsstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = native.VSStls(vsstls_input.to_native())
        vsstls.hdf_file_name = vsstls._get_hdf_file(scheme.inputs)
        vsstls._save(scheme)
        assert mock_mpi_is_root.call_count == 3
        assert os.path.isfile(vsstls.hdf_file_name)
        inspect_data = DataBase().inspect(vsstls.hdf_file_name)
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
            DataBase.ResultNames.ALPHA.value,
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(vsstls.hdf_file_name)


def test_get_free_energy_integrand(vsstls, mocker):
    arr_1d = np.ones(10)
    arr_2d = np.ones((3, 10))
    mocker.patch.object(
        DataBase,
        DataBase.read.__name__,
        return_value={
            DataBase.ResultNames.FXC_GRID.value: arr_1d,
            DataBase.ResultNames.FXC_INT.value: arr_2d,
            DataBase.ResultNames.ALPHA.value: arr_1d,
        },
    )
    fxci = vsstls.get_free_energy_integrand("dummy_file_name")
    assert np.array_equal(fxci.grid, arr_1d)
    assert np.array_equal(fxci.alpha, arr_1d)
    assert np.array_equal(fxci.integrand, arr_2d)
