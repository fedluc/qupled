import os
import pytest
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from qupled.util import HDF


@pytest.fixture
def hdf_instance():
    return HDF()


def mock_output(hdf_file_name):
    data1D = np.zeros(2)
    data2D = np.zeros((2, 2))
    pd.DataFrame(
        {
            HDF.EntryKeys.COUPLING.value: 0.0,
            HDF.EntryKeys.DEGENERACY.value: 0.0,
            HDF.EntryKeys.ERROR.value: 0.0,
            HDF.EntryKeys.THEORY.value: "theory",
            HDF.EntryKeys.RESOLUTION.value: 0.0,
            HDF.EntryKeys.CUTOFF.value: 0,
            HDF.EntryKeys.FREQUENCY_CUTOFF.value: 0,
            HDF.EntryKeys.MATSUBARA.value: 0,
        },
        index=[HDF.EntryKeys.INFO.value],
    ).to_hdf(hdf_file_name, key=HDF.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.ALPHA.value)
    pd.DataFrame(data2D).to_hdf(hdf_file_name, key=HDF.EntryKeys.ADR.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.BF.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.FXC_GRID.value)
    pd.DataFrame(data2D).to_hdf(hdf_file_name, key=HDF.EntryKeys.FXC_INT.value)
    pd.DataFrame(data2D).to_hdf(hdf_file_name, key=HDF.EntryKeys.IDR.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.RDF.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.RDF_GRID.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.SDR.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.SLFC.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.SSF.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.SSF_HF.value)
    pd.DataFrame(data1D).to_hdf(hdf_file_name, key=HDF.EntryKeys.WVG.value)


def mock_rdf_output(hdf_file_name):
    wvg_data = np.arange(0, 5, 0.1)
    ssf_data = np.ones(len(wvg_data))
    pd.DataFrame(
        {
            HDF.EntryKeys.COUPLING.value: 1.0,
            HDF.EntryKeys.DEGENERACY.value: 1.0,
            HDF.EntryKeys.ERROR.value: 0.0,
            HDF.EntryKeys.THEORY.value: "theory",
            HDF.EntryKeys.RESOLUTION.value: 0.0,
            HDF.EntryKeys.CUTOFF.value: 0,
            HDF.EntryKeys.FREQUENCY_CUTOFF.value: 0,
            HDF.EntryKeys.MATSUBARA.value: 0,
        },
        index=[HDF.EntryKeys.INFO.value],
    ).to_hdf(hdf_file_name, key=HDF.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(ssf_data).to_hdf(hdf_file_name, key=HDF.EntryKeys.SSF.value)
    pd.DataFrame(wvg_data).to_hdf(hdf_file_name, key=HDF.EntryKeys.WVG.value)


def test_entry_keys():
    expected_keys = {
        "ALPHA": "alpha",
        "ADR": "adr",
        "BF": "bf",
        "COUPLING": "coupling",
        "CUTOFF": "cutoff",
        "FREQUENCY_CUTOFF": "frequency_cutoff",
        "DEGENERACY": "degeneracy",
        "ERROR": "error",
        "FXC_GRID": "fxc_grid",
        "FXC_INT": "fxc_int",
        "MATSUBARA": "matsubara",
        "IDR": "idr",
        "INFO": "info",
        "RESOLUTION": "resolution",
        "RDF": "rdf",
        "RDF_GRID": "rdf_grid",
        "SDR": "sdr",
        "SLFC": "slfc",
        "SSF": "ssf",
        "SSF_HF": "ssf_HF",
        "THEORY": "theory",
        "WVG": "wvg",
    }
    assert len(HDF.EntryKeys) == len(expected_keys)
    for key, value in expected_keys.items():
        assert hasattr(HDF.EntryKeys, key)
        assert getattr(HDF.EntryKeys, key).value == value


def test_entry_type():
    expected_types = {
        "NUMPY": "numpy",
        "NUMPY2D": "numpy2D",
        "NUMBER": "number",
        "STRING": "string",
    }
    assert len(HDF.EntryType) == len(expected_types)
    for key, value in expected_types.items():
        assert hasattr(HDF.EntryType, key)
        assert getattr(HDF.EntryType, key).value == value


def test_set_entries(hdf_instance):
    for key, entry in hdf_instance.entries.items():
        if entry.entry_type == HDF.EntryType.NUMPY.value:
            value = np.array([1, 2, 3, 4])
        elif entry.entry_type == HDF.EntryType.NUMPY2D.value:
            value = np.array([1, 2, 3, 4]).reshape((2, 2))
        elif entry.entry_type == HDF.EntryType.NUMBER.value:
            value = 42
        elif entry.entry_type == HDF.EntryType.STRING.value:
            value = "test_value"
        else:
            assert False

        # Set value
        hdf_instance.entries[key] = value


def test_read(hdf_instance):
    hdf_file_name = "testOutput.h5"
    mock_output(hdf_file_name)
    all_hdf_entries = hdf_instance.entries.keys()
    read_data = hdf_instance.read(hdf_file_name, all_hdf_entries)
    try:
        for entry in all_hdf_entries:
            if entry in [
                HDF.EntryKeys.COUPLING.value,
                HDF.EntryKeys.DEGENERACY.value,
                HDF.EntryKeys.ERROR.value,
                HDF.EntryKeys.RESOLUTION.value,
                HDF.EntryKeys.CUTOFF.value,
                HDF.EntryKeys.FREQUENCY_CUTOFF.value,
                HDF.EntryKeys.MATSUBARA.value,
            ]:
                assert read_data[entry] == 0.0
            elif entry in [
                HDF.EntryKeys.BF.value,
                HDF.EntryKeys.FXC_GRID.value,
                HDF.EntryKeys.RDF.value,
                HDF.EntryKeys.RDF_GRID.value,
                HDF.EntryKeys.SDR.value,
                HDF.EntryKeys.SLFC.value,
                HDF.EntryKeys.SSF.value,
                HDF.EntryKeys.SSF_HF.value,
                HDF.EntryKeys.WVG.value,
                HDF.EntryKeys.ALPHA.value,
            ]:
                assert np.array_equal(read_data[entry], np.zeros(2))
            elif entry in [
                HDF.EntryKeys.ADR.value,
                HDF.EntryKeys.FXC_INT.value,
                HDF.EntryKeys.IDR.value,
            ]:
                assert np.array_equal(read_data[entry], np.zeros((2, 2)))
            elif entry == HDF.EntryKeys.THEORY.value:
                assert read_data[entry] == "theory"
            else:
                assert False
        with pytest.raises(KeyError) as excinfo:
            hdf_instance.read(hdf_file_name, ["dummyEntry"])
            assert str(excinfo.value) == "Unknown entry: dummyEntry"
    finally:
        os.remove(hdf_file_name)


def test_inspect(hdf_instance):
    hdf_file_name = "testOutput.h5"
    mock_output(hdf_file_name)
    all_hdf_entries = hdf_instance.entries.keys()
    inspect_data = hdf_instance.inspect(hdf_file_name)
    try:
        for entry in all_hdf_entries:
            assert entry in list(inspect_data.keys())
            assert inspect_data[entry] == hdf_instance.entries[entry].description
    finally:
        os.remove(hdf_file_name)


def test_plot(hdf_instance, mocker):
    hdf_file_name = "testOutput.h5"
    mock_plot_show = mocker.patch.object(plt, plt.show.__name__)
    mock_output(hdf_file_name)
    to_plot = [
        HDF.EntryKeys.RDF.value,
        HDF.EntryKeys.ADR.value,
        HDF.EntryKeys.IDR.value,
        HDF.EntryKeys.FXC_INT.value,
        HDF.EntryKeys.BF.value,
        HDF.EntryKeys.SDR.value,
        HDF.EntryKeys.SLFC.value,
        HDF.EntryKeys.SSF.value,
        HDF.EntryKeys.SSF_HF.value,
        HDF.EntryKeys.ALPHA.value,
    ]
    try:
        hdf_instance.plot(hdf_file_name, to_plot)
        assert mock_plot_show.call_count == len(to_plot)
        with pytest.raises(ValueError) as excinfo:
            hdf_instance.plot(hdf_file_name, ["dummyQuantityToPlot"])
            assert str(excinfo.value) == "Unknown quantity to plot: dummyQuantityToPlot"
    finally:
        os.remove(hdf_file_name)


def test_compute_rdf(hdf_instance):
    hdf_file_name = "testOutput.h5"
    mock_rdf_output(hdf_file_name)
    try:
        hdf_instance.compute_rdf(hdf_file_name, np.arange(0, 10, 0.1), False)
        inspect_data = hdf_instance.inspect(hdf_file_name)
        assert HDF.EntryKeys.RDF.value not in list(inspect_data.keys())
        assert HDF.EntryKeys.RDF_GRID.value not in list(inspect_data.keys())
        hdf_instance.compute_rdf(hdf_file_name, np.arange(0, 10, 0.1), True)
        inspect_data = hdf_instance.inspect(hdf_file_name)
        assert HDF.EntryKeys.RDF.value in list(inspect_data.keys())
        assert HDF.EntryKeys.RDF_GRID.value in list(inspect_data.keys())
    finally:
        os.remove(hdf_file_name)


def test_compute_internal_energy(hdf_instance):
    hdf_file_name = "testOutput.h5"
    mock_rdf_output(hdf_file_name)
    try:
        uint = hdf_instance.compute_internal_energy(hdf_file_name)
        assert uint == 0.0
    finally:
        os.remove(hdf_file_name)
