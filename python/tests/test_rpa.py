import os
import pytest
from qupled.qupled import Rpa as RpaNative
from qupled.util import Hdf
from qupled.classic import Rpa


@pytest.fixture
def rpa():
    return Rpa()


@pytest.fixture
def rpa_input():
    return Rpa.Input(1.0, 1.0)


def test_default(rpa):
    assert rpa.hdfFileName is None


def test_compute(rpa, rpa_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Rpa.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Rpa._checkStatusAndClean")
    mockSave = mocker.patch("qupled.classic.Rpa._save")
    rpa.compute(rpa_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_checkStatusAndClean(rpa, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockCheckInputs = mocker.patch("os.remove")
    rpa._checkStatusAndClean(0, "")
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 1
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        rpa._checkStatusAndClean(1, "")
    assert excinfo.value.code == "Error while solving the dielectric theory"


def test_getHdfFile(rpa, rpa_input):
    filename = rpa._getHdfFile(rpa_input)
    assert filename == "rs1.000_theta1.000_RPA.h5"


def test_save(rpa, rpa_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = RpaNative(rpa_input.toNative())
        rpa.hdfFileName = rpa._getHdfFile(scheme.inputs)
        rpa._save(scheme)
        assert mockMPIIsRoot.call_count == 1
        assert os.path.isfile(rpa.hdfFileName)
        inspectData = Hdf().inspect(rpa.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "resolution",
            "cutoff",
            "matsubara",
            "idr",
            "sdr",
            "slfc",
            "ssf",
            "ssfHF",
            "wvg",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(rpa.hdfFileName)


def test_computeRdf(rpa, mocker):
    mockMPIGetRank = mocker.patch("qupled.util.MPI.getRank", return_value=0)
    mockComputeRdf = mocker.patch("qupled.util.Hdf.computeRdf")
    rpa.computeRdf()
    assert mockMPIGetRank.call_count == 1
    assert mockComputeRdf.call_count == 1


def test_computeInternalEnergy(rpa, mocker):
    mockComputeInternalEnergy = mocker.patch("qupled.util.Hdf.computeInternalEnergy")
    rpa.computeInternalEnergy()
    assert mockComputeInternalEnergy.call_count == 1


def test_plot(rpa, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockComputeRdf = mocker.patch("qupled.classic.Rpa.computeRdf")
    mockPlot = mocker.patch("qupled.util.Hdf.plot")
    rpa.plot(["ssf", "idr"])
    assert mockMPIIsRoot.call_count == 1
    assert mockComputeRdf.call_count == 0
    assert mockPlot.call_count == 1
    rpa.plot(["ssf", "rdf"])
    assert mockMPIIsRoot.call_count == 2
    assert mockComputeRdf.call_count == 1
    assert mockPlot.call_count == 2
