import os
import pytest
from qupled.native import Stls as NativeStls
from qupled.util import Hdf, MPI
from qupled.classic import StlsIet


@pytest.fixture
def stls_iet():
    return StlsIet()


@pytest.fixture
def stls_iet_input():
    return StlsIet.Input(1.0, 1.0, "STLS-HNC")


def test_default(stls_iet):
    assert stls_iet.hdfFileName is None


def test_compute(stls_iet, stls_iet_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(StlsIet, StlsIet._compute.__name__)
    mockSave = mocker.patch.object(StlsIet, StlsIet._save.__name__)
    stls_iet.compute(stls_iet_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_iet, stls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    try:
        scheme = NativeStls(stls_iet_input.toNative())
        stls_iet.hdfFileName = stls_iet._getHdfFile(scheme.inputs)
        stls_iet._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(stls_iet.hdfFileName)
        inspectData = Hdf().inspect(stls_iet.hdfFileName)
        expectedEntries = [
            Hdf.EntryKeys.COUPLING.value,
            Hdf.EntryKeys.DEGENERACY.value,
            Hdf.EntryKeys.THEORY.value,
            Hdf.EntryKeys.ERROR.value,
            Hdf.EntryKeys.RESOLUTION.value,
            Hdf.EntryKeys.CUTOFF.value,
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value,
            Hdf.EntryKeys.MATSUBARA.value,
            Hdf.EntryKeys.BF.value,
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
        os.remove(stls_iet.hdfFileName)
