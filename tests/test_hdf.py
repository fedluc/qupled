import os
import pytest
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from qupled.util import Hdf


@pytest.fixture
def hdf_instance():
    return Hdf()


def mockOutput(hdfFileName):
    data1D = np.zeros(2)
    data2D = np.zeros((2, 2))
    pd.DataFrame(
        {
            Hdf.EntryKeys.COUPLING.value: 0.0,
            Hdf.EntryKeys.DEGENERACY.value: 0.0,
            Hdf.EntryKeys.ERROR.value: 0.0,
            Hdf.EntryKeys.THEORY.value: "theory",
            Hdf.EntryKeys.RESOLUTION.value: 0.0,
            Hdf.EntryKeys.CUTOFF.value: 0,
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value: 0,
            Hdf.EntryKeys.MATSUBARA.value: 0,
        },
        index=[Hdf.EntryKeys.INFO.value],
    ).to_hdf(hdfFileName, key=Hdf.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.ALPHA.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=Hdf.EntryKeys.ADR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.BF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.FXC_GRID.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=Hdf.EntryKeys.FXCI.value)
    pd.DataFrame(data2D).to_hdf(hdfFileName, key=Hdf.EntryKeys.IDR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.RDF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.RDF_GRID.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.SDR.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.SLFC.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.SSF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.SSF_HF.value)
    pd.DataFrame(data1D).to_hdf(hdfFileName, key=Hdf.EntryKeys.WVG.value)


def mockRdfOutput(hdfFileName):
    wvgData = np.arange(0, 5, 0.1)
    ssfData = np.ones(len(wvgData))
    pd.DataFrame(
        {
            Hdf.EntryKeys.COUPLING.value: 1.0,
            Hdf.EntryKeys.DEGENERACY.value: 1.0,
            Hdf.EntryKeys.ERROR.value: 0.0,
            Hdf.EntryKeys.THEORY.value: "theory",
            Hdf.EntryKeys.RESOLUTION.value: 0.0,
            Hdf.EntryKeys.CUTOFF.value: 0,
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value: 0,
            Hdf.EntryKeys.MATSUBARA.value: 0,
        },
        index=[Hdf.EntryKeys.INFO.value],
    ).to_hdf(hdfFileName, key=Hdf.EntryKeys.INFO.value, mode="w")
    pd.DataFrame(ssfData).to_hdf(hdfFileName, key=Hdf.EntryKeys.SSF.value)
    pd.DataFrame(wvgData).to_hdf(hdfFileName, key=Hdf.EntryKeys.WVG.value)


def test_entry_keys():
    expected_keys = {
        "ALPHA": "alpha",
        "ADR": "adr",
        "BF": "bf",
        "COUPLING": "coupling",
        "CUTOFF": "cutoff",
        "FREQUENCY_CUTOFF": "frequencyCutoff",
        "DEGENERACY": "degeneracy",
        "ERROR": "error",
        "FXC_GRID": "fxcGrid",
        "FXCI": "fxci",
        "MATSUBARA": "matsubara",
        "IDR": "idr",
        "INFO": "info",
        "RESOLUTION": "resolution",
        "RDF": "rdf",
        "RDF_GRID": "rdfGrid",
        "SDR": "sdr",
        "SLFC": "slfc",
        "SSF": "ssf",
        "SSF_HF": "ssfHF",
        "THEORY": "theory",
        "WVG": "wvg",
    }
    assert len(Hdf.EntryKeys) == len(expected_keys)
    for key, value in expected_keys.items():
        assert hasattr(Hdf.EntryKeys, key)
        assert getattr(Hdf.EntryKeys, key).value == value


def test_entry_type():
    expected_types = {
        "NUMPY": "numpy",
        "NUMPY2D": "numpy2D",
        "NUMBER": "number",
        "STRING": "string",
    }
    assert len(Hdf.EntryType) == len(expected_types)
    for key, value in expected_types.items():
        assert hasattr(Hdf.EntryType, key)
        assert getattr(Hdf.EntryType, key).value == value


def test_set_entries(hdf_instance):
    for key, entry in hdf_instance.entries.items():
        if entry.entryType == Hdf.EntryType.NUMPY.value:
            value = np.array([1, 2, 3, 4])
        elif entry.entryType == Hdf.EntryType.NUMPY2D.value:
            value = np.array([1, 2, 3, 4]).reshape((2, 2))
        elif entry.entryType == Hdf.EntryType.NUMBER.value:
            value = 42
        elif entry.entryType == Hdf.EntryType.STRING.value:
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
                Hdf.EntryKeys.COUPLING.value,
                Hdf.EntryKeys.DEGENERACY.value,
                Hdf.EntryKeys.ERROR.value,
                Hdf.EntryKeys.RESOLUTION.value,
                Hdf.EntryKeys.CUTOFF.value,
                Hdf.EntryKeys.FREQUENCY_CUTOFF.value,
                Hdf.EntryKeys.MATSUBARA.value,
            ]:
                assert readData[entry] == 0.0
            elif entry in [
                Hdf.EntryKeys.BF.value,
                Hdf.EntryKeys.FXC_GRID.value,
                Hdf.EntryKeys.RDF.value,
                Hdf.EntryKeys.RDF_GRID.value,
                Hdf.EntryKeys.SDR.value,
                Hdf.EntryKeys.SLFC.value,
                Hdf.EntryKeys.SSF.value,
                Hdf.EntryKeys.SSF_HF.value,
                Hdf.EntryKeys.WVG.value,
                Hdf.EntryKeys.ALPHA.value,
            ]:
                assert np.array_equal(readData[entry], np.zeros(2))
            elif entry in [
                Hdf.EntryKeys.ADR.value,
                Hdf.EntryKeys.FXCI.value,
                Hdf.EntryKeys.IDR.value,
            ]:
                assert np.array_equal(readData[entry], np.zeros((2, 2)))
            elif entry == Hdf.EntryKeys.THEORY.value:
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
        Hdf.EntryKeys.RDF.value,
        Hdf.EntryKeys.ADR.value,
        Hdf.EntryKeys.IDR.value,
        Hdf.EntryKeys.FXCI.value,
        Hdf.EntryKeys.BF.value,
        Hdf.EntryKeys.SDR.value,
        Hdf.EntryKeys.SLFC.value,
        Hdf.EntryKeys.SSF.value,
        Hdf.EntryKeys.SSF_HF.value,
        Hdf.EntryKeys.ALPHA.value,
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
        hdf_instance.computeRdf(hdfFileName, np.arange(0, 10, 0.1), False)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert Hdf.EntryKeys.RDF.value not in list(inspectData.keys())
        assert Hdf.EntryKeys.RDF_GRID.value not in list(inspectData.keys())
        hdf_instance.computeRdf(hdfFileName, np.arange(0, 10, 0.1), True)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert Hdf.EntryKeys.RDF.value in list(inspectData.keys())
        assert Hdf.EntryKeys.RDF_GRID.value in list(inspectData.keys())
    finally:
        os.remove(hdfFileName)


def test_computeInternalEnergy(hdf_instance):
    hdfFileName = "testOutput.h5"
    mockRdfOutput(hdfFileName)
    try:
        uint = hdf_instance.computeInternalEnergy(hdfFileName)
        assert uint == 0.0
    finally:
        os.remove(hdfFileName)
