import os
import pytest
import numpy as np
from qupled.native import Qstls as NativeQstls
from qupled.util import Hdf, MPI
from qupled.quantum import Qstls


@pytest.fixture
def qstls():
    return Qstls()


@pytest.fixture
def qstls_input():
    return Qstls.Input(1.0, 1.0)


def test_default(qstls):
    assert qstls.hdfFileName is None


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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    try:
        scheme = NativeQstls(qstls_input.toNative())
        qstls.hdfFileName = qstls._getHdfFile(scheme.inputs)
        qstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(qstls.hdfFileName)
        inspectData = Hdf().inspect(qstls.hdfFileName)
        expectedEntries = [
            Hdf.EntryKeys.COUPLING.value,
            Hdf.EntryKeys.DEGENERACY.value,
            Hdf.EntryKeys.THEORY.value,
            Hdf.EntryKeys.ERROR.value,
            Hdf.EntryKeys.RESOLUTION.value,
            Hdf.EntryKeys.CUTOFF.value,
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value,
            Hdf.EntryKeys.MATSUBARA.value,
            Hdf.EntryKeys.ADR.value,
            Hdf.EntryKeys.IDR.value,
            Hdf.EntryKeys.SDR.value,
            Hdf.EntryKeys.SLFC.value,
            Hdf.EntryKeys.SSF.value,
            Hdf.EntryKeys.SSF_HF.value,
            Hdf.EntryKeys.WVG.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qstls.hdfFileName)


def test_getInitialGuess(mocker):
    arr = np.ones(10)
    mocker.patch.object(
        Hdf,
        Hdf.read.__name__,
        return_value={
            Hdf.EntryKeys.WVG.value: arr,
            Hdf.EntryKeys.SSF.value: arr,
            Hdf.EntryKeys.ADR.value: arr,
            Hdf.EntryKeys.MATSUBARA.value: 10,
        },
    )
    guess = Qstls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.ssf, arr)
    assert np.array_equal(guess.adr, arr)
    assert np.array_equal(guess.matsubara, 10)
