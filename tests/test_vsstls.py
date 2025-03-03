import os
import pytest
import numpy as np
from qupled import native
from qupled.util import HDF, MPI
from qupled.classic import VSStls


@pytest.fixture
def vsstls():
    return VSStls()


@pytest.fixture
def vsstls_input():
    return VSStls.Input(1.0, 1.0)


def test_default(vsstls):
    assert vsstls.hdfFileName is None


def test_compute(vsstls, vsstls_input, mocker):
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(VSStls, VSStls._compute.__name__)
    mockSave = mocker.patch.object(VSStls, VSStls._save.__name__)
    vsstls.compute(vsstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1


def test_save(vsstls, vsstls_input, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = native.VSStls(vsstls_input.toNative())
        vsstls.hdfFileName = vsstls._getHdfFile(scheme.inputs)
        vsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(vsstls.hdfFileName)
        inspectData = HDF().inspect(vsstls.hdfFileName)
        expectedEntries = [
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
            HDF.EntryKeys.FXC_GRID.value,
            HDF.EntryKeys.FXCI.value,
            HDF.EntryKeys.ALPHA.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(vsstls.hdfFileName)


def test_getFreeEnergyIntegrand(vsstls, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mocker.patch.object(
        HDF,
        HDF.read.__name__,
        return_value={
            HDF.EntryKeys.FXC_GRID.value: arr1D,
            HDF.EntryKeys.FXCI.value: arr2D,
            HDF.EntryKeys.ALPHA.value: arr1D,
        },
    )
    fxci = vsstls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxci.grid, arr1D)
    assert np.array_equal(fxci.alpha, arr1D)
    assert np.array_equal(fxci.integrand, arr2D)
