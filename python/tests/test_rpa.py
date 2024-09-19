import os
import pytest
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import Rpa, ClassicScheme


@pytest.fixture
def rpa_instance():
    inputs = qp.RpaInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "RPA"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    return Rpa(Rpa.Input(1.0, 1.0))


def test_default(rpa_instance):
    issubclass(Rpa, ClassicScheme)
    issubclass(Rpa, qp.Rpa)
    assert all(x == y for x, y in zip(rpa_instance.allowedTheories, ["RPA"]))
    assert rpa_instance.hdfFileName is None


def test_compute(rpa_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.Rpa._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Rpa.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Rpa._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.Rpa._setHdfFile")
    mockSave = mocker.patch("qupled.classic.Rpa._save")
    rpa_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_checkStatusAndClean(rpa_instance, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockCheckInputs = mocker.patch("os.remove")
    rpa_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert mockMPIIsRoot.call_count == 1
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        rpa_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"


def test_setHdfFile(rpa_instance):
    rpa_instance._setHdfFile()
    assert rpa_instance.hdfFileName == "rs1.000_theta1.000_RPA.h5"


def test_save(rpa_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    rpa_instance._setHdfFile()
    try:
        rpa_instance._save()
        assert mockMPIIsRoot.call_count == 1
        assert os.path.isfile(rpa_instance.hdfFileName)
        inspectData = Hdf().inspect(rpa_instance.hdfFileName)
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
        os.remove(rpa_instance.hdfFileName)


def test_computeRdf(rpa_instance, mocker):
    mockMPIGetRank = mocker.patch("qupled.util.MPI.getRank", return_value=0)
    mockComputeRdf = mocker.patch("qupled.util.Hdf.computeRdf")
    rpa_instance.computeRdf()
    assert mockMPIGetRank.call_count == 1
    assert mockComputeRdf.call_count == 1


def test_computeInternalEnergy(rpa_instance, mocker):
    mockComputeInternalEnergy = mocker.patch("qupled.qupled.computeInternalEnergy")
    rpa_instance.computeInternalEnergy()
    assert mockComputeInternalEnergy.call_count == 1


def test_plot(rpa_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockComputeRdf = mocker.patch("qupled.classic.Rpa.computeRdf")
    mockPlot = mocker.patch("qupled.util.Hdf.plot")
    rpa_instance.plot(["ssf", "idr"])
    assert mockMPIIsRoot.call_count == 1
    assert mockComputeRdf.call_count == 0
    assert mockPlot.call_count == 1
    rpa_instance.plot(["ssf", "rdf"])
    assert mockMPIIsRoot.call_count == 2
    assert mockComputeRdf.call_count == 1
    assert mockPlot.call_count == 2
