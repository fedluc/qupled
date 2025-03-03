import os
import pytest
import numpy as np
from qupled.native import Qstls as NativeQstls
from qupled.util import HDF, MPI
from qupled.quantum import Qstls


@pytest.fixture
def qstls():
    return Qstls()


@pytest.fixture
def qstls_input():
    return Qstls.Input(1.0, 1.0)


def test_default(qstls):
    assert qstls.hdf_file_name is None


def test_compute(qstls, qstls_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(Qstls, Qstls._compute.__name__)
    mockSave = mocker.patch.object(Qstls, Qstls._save.__name__)
    qstls.compute(qstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1


def test_save(qstls, qstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeQstls(qstls_input.toNative())
        qstls.hdf_file_name = qstls._get_hdf_file(scheme.inputs)
        qstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(qstls.hdf_file_name)
        inspectData = HDF().inspect(qstls.hdf_file_name)
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
            HDF.EntryKeys.SSF.value,
            HDF.EntryKeys.SSF_HF.value,
            HDF.EntryKeys.WVG.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qstls.hdf_file_name)


def test_getInitialGuess(mocker):
    arr = np.ones(10)
    mocker.patch.object(
        HDF,
        HDF.read.__name__,
        return_value={
            HDF.EntryKeys.WVG.value: arr,
            HDF.EntryKeys.SSF.value: arr,
            HDF.EntryKeys.ADR.value: arr,
            HDF.EntryKeys.MATSUBARA.value: 10,
        },
    )
    guess = Qstls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.ssf, arr)
    assert np.array_equal(guess.adr, arr)
    assert np.array_equal(guess.matsubara, 10)
