import os
import pytest
from qupled.native import Rpa as NativeRpa
from qupled.util import HDF, MPI
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
    mockMPITime = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mockMPIBarrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mockCompute = mocker.patch.object(Rpa, Rpa._compute.__name__)
    mockSave = mocker.patch.object(Rpa, Rpa._save.__name__)
    rpa.compute(rpa_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1


def test_checkStatusAndClean(rpa, mocker, capsys):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mocker.patch.object(os, os.remove.__name__)
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
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeRpa(rpa_input.toNative())
        rpa.hdfFileName = rpa._getHdfFile(scheme.inputs)
        rpa._save(scheme)
        assert mockMPIIsRoot.call_count == 1
        assert os.path.isfile(rpa.hdfFileName)
        inspectData = HDF().inspect(rpa.hdfFileName)
        expectedEntries = [
            HDF.EntryKeys.COUPLING.value,
            HDF.EntryKeys.DEGENERACY.value,
            HDF.EntryKeys.THEORY.value,
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
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(rpa.hdfFileName)


def test_computeRdf(rpa, mocker):
    mockMPIGetRank = mocker.patch.object(MPI, MPI.rank.__name__, return_value=0)
    mockComputeRdf = mocker.patch.object(HDF, HDF.compute_rdf.__name__)
    rpa.computeRdf()
    assert mockMPIGetRank.call_count == 1
    assert mockComputeRdf.call_count == 1


def test_computeInternalEnergy(rpa, mocker):
    mockComputeInternalEnergy = mocker.patch.object(
        HDF, HDF.compute_internal_energy.__name__
    )
    rpa.computeInternalEnergy()
    assert mockComputeInternalEnergy.call_count == 1


def test_plot(rpa, mocker):
    mockMPIIsRoot = mocker.patch.object(MPI, MPI.is_root.__name__)
    mockComputeRdf = mocker.patch.object(HDF, HDF.compute_rdf.__name__)
    mockPlot = mocker.patch.object(HDF, HDF.plot.__name__)
    rpa.plot([HDF.EntryKeys.SSF.value, HDF.EntryKeys.IDR.value])
    assert mockMPIIsRoot.call_count == 1
    assert mockComputeRdf.call_count == 0
    assert mockPlot.call_count == 1
    rpa.plot([HDF.EntryKeys.SSF.value, HDF.EntryKeys.RDF.value])
    assert mockMPIIsRoot.call_count == 2
    assert mockComputeRdf.call_count == 1
    assert mockPlot.call_count == 2
