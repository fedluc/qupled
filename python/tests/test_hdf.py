import os
import pytest
import numpy as np
import pandas as pd
from qupled.util import Hdf


@pytest.fixture
def hdf_instance():
    return Hdf()


def mockOutput(hdfFileName):
    data1D = np.zeros(2)
    data2D = np.zeros((2, 2))
    pd.DataFrame(
        {
            "coupling": 0.0,
            "degeneracy": 0.0,
            "error": 0.0,
            "theory": "theory",
            "resolution": 0.0,
            "cutoff": 0,
            "matsubara": 0,
        },
        index=["info"],
    ).to_hdf(hdfFileName, key="info", mode="w")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="alpha")
    pd.DataFrame(data2D).to_hdf(hdfFileName, key="adr")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="bf")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="fxcGrid")
    pd.DataFrame(data2D).to_hdf(hdfFileName, key="fxci")
    pd.DataFrame(data2D).to_hdf(hdfFileName, key="idr")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="rdf")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="rdfGrid")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="sdr")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="slfc")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="ssf")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="ssfHF")
    pd.DataFrame(data1D).to_hdf(hdfFileName, key="wvg")


def mockRdfOutput(hdfFileName):
    wvgData = np.arange(0, 5, 0.1)
    ssfData = np.ones(len(wvgData))
    pd.DataFrame(
        {
            "coupling": 1.0,
            "degeneracy": 1.0,
            "error": 0.0,
            "theory": "theory",
            "resolution": 0.0,
            "cutoff": 0,
            "matsubara": 0,
        },
        index=["info"],
    ).to_hdf(hdfFileName, key="info", mode="w")
    pd.DataFrame(ssfData).to_hdf(hdfFileName, key="ssf")
    pd.DataFrame(wvgData).to_hdf(hdfFileName, key="wvg")


def test_set_entries(hdf_instance):
    for key, entry in hdf_instance.entries.items():
        if entry.entryType == "numpy":
            value = np.array([1, 2, 3, 4])
        elif entry.entryType == "numpy2D":
            value = np.array([1, 2, 3, 4]).reshape((2, 2))
        elif entry.entryType == "number":
            value = 42
        elif entry.entryType == "string":
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
                "coupling",
                "degeneracy",
                "error",
                "resolution",
                "cutoff",
                "matsubara",
            ]:
                assert readData[entry] == 0.0
            elif entry in [
                "bf",
                "fxcGrid",
                "rdf",
                "rdfGrid",
                "sdr",
                "slfc",
                "ssf",
                "ssfHF",
                "wvg",
                "alpha",
            ]:
                assert np.array_equal(readData[entry], np.zeros(2))
            elif entry in ["adr", "fxci", "idr"]:
                assert np.array_equal(readData[entry], np.zeros((2, 2)))
            elif entry in ["theory"]:
                assert readData[entry] == "theory"
            else:
                assert False
        with pytest.raises(SystemExit) as excinfo:
            hdf_instance.read(hdfFileName, "dummyEntry")
        assert excinfo.value.code == "Unknown entry"
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
    mockPlotShow = mocker.patch("matplotlib.pyplot.show")
    mockOutput(hdfFileName)
    toPlot = ["rdf", "adr", "idr", "fxci", "bf", "sdr", "slfc", "ssf", "ssfHF", "alpha"]
    try:
        hdf_instance.plot(hdfFileName, toPlot)
        assert mockPlotShow.call_count == len(toPlot)
        with pytest.raises(SystemExit) as excinfo:
            hdf_instance.plot(hdfFileName, "dummyQuantityToPlot")
        assert excinfo.value.code == "Unknown quantity to plot"
    finally:
        os.remove(hdfFileName)


def test_computeRdf(hdf_instance):
    hdfFileName = "testOutput.h5"
    mockRdfOutput(hdfFileName)
    try:
        hdf_instance.computeRdf(hdfFileName, np.arange(0, 10, 0.1), False)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert "rdf" not in list(inspectData.keys())
        assert "rdfGrid" not in list(inspectData.keys())
        hdf_instance.computeRdf(hdfFileName, np.arange(0, 10, 0.1), True)
        inspectData = hdf_instance.inspect(hdfFileName)
        assert "rdf" in list(inspectData.keys())
        assert "rdfGrid" in list(inspectData.keys())
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
