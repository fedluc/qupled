import os
import pytest
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from qupled.util import HDF


@pytest.fixture
def hdf_instance():
    return HDF()


def mockOutput(hdfFileName):
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
    ).to_hdf(hdfFileName, key=HDF.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.ALPHA.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=HDF.EntryKeys.ADR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.BF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.FXC_GRID.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=HDF.EntryKeys.FXCI.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=HDF.EntryKeys.IDR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.RDF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.RDF_GRID.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.SDR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.SLFC.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.SSF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.SSF_HF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=HDF.EntryKeys.WVG.value)


def mockRdfOutput(hdfFileName):
    wvgData = np.arange(0, 5, 0.1)
    ssfData = np.ones(len(wvgData))
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
    ).to_hdf(hdfFileName, key=HDF.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(ssfData).to_hdf(hdfFileName, key=HDF.EntryKeys.SSF.value)
    pd.DataFrame(wvgData).to_hdf(hdfFileName, key=HDF.EntryKeys.WVG.value)


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
        "FXCI": "fxc_int",
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
    hdfFileName = "testOutput.h5"
    mockOutput(hdfFileName)
    allHdfEntries = hdf_instance.entries.keys()
    readData = hdf_instance.read(hdfFileName, allHdfEntries)
    try:
        for entry in allHdfEntries:
            if entry in [
                HDF.EntryKeys.COUPLING.value,
                HDF.EntryKeys.DEGENERACY.value,
                HDF.EntryKeys.ERROR.value,
                HDF.EntryKeys.RESOLUTION.value,
                HDF.EntryKeys.CUTOFF.value,
                HDF.EntryKeys.FREQUENCY_CUTOFF.value,
                HDF.EntryKeys.MATSUBARA.value,
            ]:
                assert readData[entry] == 0.0
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
                assert np.array_equal(readData[entry], np.zeros(2))
            elif entry in [
                HDF.EntryKeys.ADR.value,
                HDF.EntryKeys.FXCI.value,
                HDF.EntryKeys.IDR.value,
            ]:
                assert np.array_equal(readData[entry], np.zeros((2, 2)))
            elif entry == HDF.EntryKeys.THEORY.value:
                assert readData[entry] == "theory"
            else:
                assert False
        with pytest.raises(KeyError) as excinfo:
            hdf_instance.read(hdfFileName, ["dummyEntry"])
            assert str(excinfo.value) == "Unknown entry: dummyEntry"
    finally:
        os.remove(hdfFileName)


def test_inspect(hdf_instance):
    hdfFileName = "testOutput.h5"
    mockOutput(hdfFileName)
    allHdfEntries = hdf_instance.entries.keys()
    inspectData = hdf_instance.inspect(hdfFileName)
    try:
        for entry in allHdfEntries:
            assert entry in list(inspectData.keys())
            assert inspectData[entry] == hdf_instance.entries[entry].description
    finally:
        os.remove(hdfFileName)


def test_plot(hdf_instance, mocker):
    hdfFileName = "testOutput.h5"
    mockPlotShow = mocker.patch.object(plt, plt.show.__name__)
    mockOutput(hdfFileName)
    toPlot = [
        HDF.EntryKeys.RDF.value,
        HDF.EntryKeys.ADR.value,
        HDF.EntryKeys.IDR.value,
        HDF.EntryKeys.FXCI.value,
        HDF.EntryKeys.BF.value,
        HDF.EntryKeys.SDR.value,
        HDF.EntryKeys.SLFC.value,
        HDF.EntryKeys.SSF.value,
        HDF.EntryKeys.SSF_HF.value,
        HDF.EntryKeys.ALPHA.value,
    ]
    try:
        hdf_instance.plot(hdfFileName, toPlot)
        assert mockPlotShow.call_count == len(toPlot)
        with pytest.raises(ValueError) as excinfo:
            hdf_instance.plot(hdfFileName, ["dummyQuantityToPlot"])
            assert str(excinfo.value) == "Unknown quantity to plot: dummyQuantityToPlot"
    finally:
        os.remove(hdfFileName)


def test_computeRdf(hdf_instance):
    hdfFileName = "testOutput.h5"
    mockRdfOutput(hdfFileName)
    try:
        hdf_instance.compute_rdf(hdfFileName, np.arange(0, 10, 0.1), False)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert HDF.EntryKeys.RDF.value not in list(inspectData.keys())
        assert HDF.EntryKeys.RDF_GRID.value not in list(inspectData.keys())
        hdf_instance.compute_rdf(hdfFileName, np.arange(0, 10, 0.1), True)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert HDF.EntryKeys.RDF.value in list(inspectData.keys())
        assert HDF.EntryKeys.RDF_GRID.value in list(inspectData.keys())
    finally:
        os.remove(hdfFileName)


def test_computeInternalEnergy(hdf_instance):
    hdfFileName = "testOutput.h5"
    mockRdfOutput(hdfFileName)
    try:
        uint = hdf_instance.compute_internal_energy(hdfFileName)
        assert uint == 0.0
    finally:
        os.remove(hdfFileName)
