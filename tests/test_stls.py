import os
import pytest
import numpy as np
from qupled.native import Stls as NativeStls
from qupled.util import Hdf, MPI
from qupled.classic import Stls


@pytest.fixture
def stls():
    return Stls()


@pytest.fixture
def stls_input():
    return Stls.Input(1.0, 1.0)


def test_default(stls):
    assert stls.hdfFileName is None


def test_compute(stls, stls_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(Stls, Stls._compute.__name__)
    mockSave = mocker.patch.object(Stls, Stls._save.__name__)
    stls.compute(stls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls, stls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    try:
        scheme = NativeStls(stls_input.toNative())
        stls.hdfFileName = stls._getHdfFile(scheme.inputs)
        stls._save(scheme)
        assert mockMPIIsRoot.call_count == 2
        assert os.path.isfile(stls.hdfFileName)
        inspectData = Hdf().inspect(stls.hdfFileName)
        expectedEntries = [
            Hdf.EntryKeys.COUPLING.value,
            Hdf.EntryKeys.DEGENERACY.value,
            Hdf.EntryKeys.THEORY.value,
            Hdf.EntryKeys.ERROR.value,
            Hdf.EntryKeys.RESOLUTION.value,
            Hdf.EntryKeys.CUTOFF.value,
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value,
            Hdf.EntryKeys.MATSUBARA.value,
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
        os.remove(stls.hdfFileName)


def test_getInitialGuess(mocker):
    arr = np.ones(10)
    mocker.patch.object(
        Hdf,
        Hdf.read.__name__,
        return_value={Hdf.EntryKeys.WVG.value: arr, Hdf.EntryKeys.SLFC.value: arr},
    )
    guess = Stls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.slfc, arr)
