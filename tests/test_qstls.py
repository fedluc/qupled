import os
import pytest
import numpy as np
from qupled.native import Qstls as NativeQstls
from qupled.util import DataBase, MPI
from qupled.qstls import Qstls


@pytest.fixture
def qstls():
    return Qstls()


@pytest.fixture
def qstls_input():
    return Qstls.Input(1.0, 1.0)


def test_default(qstls):
    assert qstls.hdf_file_name is None


def test_compute(qstls, qstls_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(Qstls, Qstls._compute.__name__)
    mock_save = mocker.patch.object(Qstls, Qstls._save.__name__)
    qstls.compute(qstls_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1


def test_save(qstls, qstls_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQstls(qstls_input.to_native())
        qstls.hdf_file_name = qstls._get_hdf_file(scheme.inputs)
        qstls._save(scheme)
        assert mock_mpi_is_root.call_count == 3
        assert os.path.isfile(qstls.hdf_file_name)
        inspect_data = DataBase().inspect(qstls.hdf_file_name)
        expected_entries = [
            DataBase.ResultNames.COUPLING.value,
            DataBase.ResultNames.DEGENERACY.value,
            DataBase.ResultNames.THEORY.value,
            DataBase.ResultNames.ERROR.value,
            DataBase.ResultNames.RESOLUTION.value,
            DataBase.ResultNames.CUTOFF.value,
            DataBase.ResultNames.FREQUENCY_CUTOFF.value,
            DataBase.ResultNames.MATSUBARA.value,
            DataBase.ResultNames.ADR.value,
            DataBase.ResultNames.IDR.value,
            DataBase.ResultNames.SDR.value,
            DataBase.ResultNames.SLFC.value,
            DataBase.ResultNames.SSF.value,
            DataBase.ResultNames.SSF_HF.value,
            DataBase.ResultNames.WVG.value,
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(qstls.hdf_file_name)


def test_get_initial_guess(mocker):
    arr = np.ones(10)
    mocker.patch.object(
        DataBase,
        DataBase.read.__name__,
        return_value={
            DataBase.ResultNames.WVG.value: arr,
            DataBase.ResultNames.SSF.value: arr,
            DataBase.ResultNames.ADR.value: arr,
            DataBase.ResultNames.MATSUBARA.value: 10,
        },
    )
    guess = Qstls.get_initial_guess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.ssf, arr)
    assert np.array_equal(guess.adr, arr)
    assert np.array_equal(guess.matsubara, 10)
