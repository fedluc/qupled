import os
import pytest
import numpy as np
from qupled.qupled import VSStls as VSStlsNative
from qupled.util import Hdf
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
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.VSStls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.VSStls._checkStatusAndClean")
    mockSave = mocker.patch("qupled.classic.VSStls._save")
    vsstls.compute(vsstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(vsstls, vsstls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = VSStlsNative(vsstls_input.toNative())
        vsstls.hdfFileName = vsstls._getHdfFile(scheme.inputs)
        vsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(vsstls.hdfFileName)
        inspectData = Hdf().inspect(vsstls.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "idr",
            "sdr",
            "slfc",
            "ssf",
            "ssfHF",
            "wvg",
            "fxcGrid",
            "fxci",
            "alpha",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(vsstls.hdfFileName)


def test_getFreeEnergyIntegrand(vsstls, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"fxcGrid": arr1D, "fxci": arr2D, "alpha": arr1D},
    )
    fxci = vsstls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxci.grid, arr1D)
    assert np.array_equal(fxci.alpha, arr1D)
    assert np.array_equal(fxci.integrand, arr2D)
