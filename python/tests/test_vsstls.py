import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import VSStls, ClassicScheme


@pytest.fixture
def vsstls_instance():
    return VSStls(VSStls.Input(1.0, 1.0))


def test_default(vsstls_instance):
    assert issubclass(VSStls, ClassicScheme)
    assert issubclass(VSStls, qp.VSStls)
    assert all(x == y for x, y in zip(vsstls_instance.allowedTheories, ["VSSTLS"]))
    assert vsstls_instance.hdfFileName is None


def test_compute(vsstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.VSStls._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.VSStls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.VSStls._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.VSStls._setHdfFile")
    mockSave = mocker.patch("qupled.classic.VSStls._save")
    vsstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_save(vsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    vsstls_instance._setHdfFile()
    try:
        vsstls_instance._save()
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(vsstls_instance.hdfFileName)
        inspectData = Hdf().inspect(vsstls_instance.hdfFileName)
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
        os.remove(vsstls_instance.hdfFileName)


def test_getFreeEnergyIntegrand(vsstls_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"fxcGrid": arr1D, "fxci": arr2D, "alpha": arr1D},
    )
    fxci = vsstls_instance.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxci.grid, arr1D)
    assert np.array_equal(fxci.alpha, arr1D)
    assert np.array_equal(fxci.integrand, arr2D)
