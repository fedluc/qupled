import os

import pytest

from qupled.rpa import Rpa
from qupled.native import Rpa as NativeRpa
from qupled.util import HDF, MPI


@pytest.fixture
def rpa():
    return Rpa()


@pytest.fixture
def rpa_input():
    return Rpa.Input(1.0, 1.0)


def test_default(rpa):
    assert rpa.hdf_file_name is None


def test_compute(rpa, rpa_input, mocker):
    mock_mpi_time = mocker.patch.object(MPI, MPI.timer.__name__, return_value=0)
    mock_mpi_barrier = mocker.patch.object(MPI, MPI.barrier.__name__)
    mock_compute = mocker.patch.object(Rpa, Rpa._compute.__name__)
    mock_save = mocker.patch.object(Rpa, Rpa._save.__name__)
    rpa.compute(rpa_input)
    assert mock_mpi_time.call_count == 2
    assert mock_mpi_barrier.call_count == 1
    assert mock_compute.call_count == 1
    assert mock_save.call_count == 1


def test_check_status_and_clean(rpa, mocker, capsys):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mocker.patch.object(os, os.remove.__name__)
    rpa._check_status_and_clean(0, "")
    captured = capsys.readouterr()
    assert mock_mpi_is_root.call_count == 1
    assert "Dielectric theory solved successfully!\n" in captured.out
    with pytest.raises(SystemExit) as excinfo:
        rpa._check_status_and_clean(1, "")
    assert excinfo.value.code == "Error while solving the dielectric theory"


def test_get_hdf_file(rpa, rpa_input):
    filename = rpa._get_hdf_file(rpa_input)
    assert filename == "rs1.000_theta1.000_RPA.h5"


def test_save(rpa, rpa_input, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    try:
        scheme = NativeRpa(rpa_input.to_native())
        rpa.hdf_file_name = rpa._get_hdf_file(scheme.inputs)
        rpa._save(scheme)
        assert mock_mpi_is_root.call_count == 1
        assert os.path.isfile(rpa.hdf_file_name)
        inspect_data = HDF().inspect(rpa.hdf_file_name)
        expected_entries = [
            HDF.ResultNames.COUPLING.value,
            HDF.ResultNames.DEGENERACY.value,
            HDF.ResultNames.THEORY.value,
            HDF.ResultNames.RESOLUTION.value,
            HDF.ResultNames.CUTOFF.value,
            HDF.ResultNames.FREQUENCY_CUTOFF.value,
            HDF.ResultNames.MATSUBARA.value,
            HDF.ResultNames.IDR.value,
            HDF.ResultNames.SDR.value,
            HDF.ResultNames.SLFC.value,
            HDF.ResultNames.SSF.value,
            HDF.ResultNames.SSF_HF.value,
            HDF.ResultNames.WVG.value,
        ]
        for entry in expected_entries:
            assert entry in inspect_data
    finally:
        os.remove(rpa.hdf_file_name)


def test_compute_rdf(rpa, mocker):
    mock_mpi_get_rank = mocker.patch.object(MPI, MPI.rank.__name__, return_value=0)
    mock_compute_rdf = mocker.patch.object(HDF, HDF.compute_rdf.__name__)
    rpa.compute_rdf()
    assert mock_mpi_get_rank.call_count == 1
    assert mock_compute_rdf.call_count == 1


def test_compute_internal_energy(rpa, mocker):
    mock_compute_internal_energy = mocker.patch.object(
        HDF, HDF.compute_internal_energy.__name__
    )
    rpa.compute_internal_energy()
    assert mock_compute_internal_energy.call_count == 1


def test_plot(rpa, mocker):
    mock_mpi_is_root = mocker.patch.object(MPI, MPI.is_root.__name__)
    mock_compute_rdf = mocker.patch.object(HDF, HDF.compute_rdf.__name__)
    mock_plot = mocker.patch.object(HDF, HDF.plot.__name__)
    rpa.plot([HDF.ResultNames.SSF.value, HDF.ResultNames.IDR.value])
    assert mock_mpi_is_root.call_count == 1
    assert mock_compute_rdf.call_count == 0
    assert mock_plot.call_count == 1
    rpa.plot([HDF.ResultNames.SSF.value, HDF.ResultNames.RDF.value])
    assert mock_mpi_is_root.call_count == 2
    assert mock_compute_rdf.call_count == 1
    assert mock_plot.call_count == 2
