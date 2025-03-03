import os
import pytest
import numpy as np
from qupled import native
from qupled.util import Hdf, MPI
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.isRoot.__name__)
    try:
        scheme = native.VSStls(vsstls_input.toNative())
        vsstls.hdfFileName = vsstls._getHdfFile(scheme.inputs)
        vsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(vsstls.hdfFileName)
        inspectData = Hdf().inspect(vsstls.hdfFileName)
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
            Hdf.EntryKeys.FXC_GRID.value,
            Hdf.EntryKeys.FXCI.value,
            Hdf.EntryKeys.ALPHA.value,
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(vsstls.hdfFileName)


def test_getFreeEnergyIntegrand(vsstls, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mocker.patch.object(
        Hdf,
        Hdf.read.__name__,
        return_value={
            Hdf.EntryKeys.FXC_GRID.value: arr1D,
            Hdf.EntryKeys.FXCI.value: arr2D,
            Hdf.EntryKeys.ALPHA.value: arr1D,
        },
    )
    fxci = vsstls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxci.grid, arr1D)
    assert np.array_equal(fxci.alpha, arr1D)
    assert np.array_equal(fxci.integrand, arr2D)
