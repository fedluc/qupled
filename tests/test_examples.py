import os
import glob
import importlib
import pytest


def cleanExampleFiles():
    for fileExtension in ["h5", "bin", "zip"]:
        fileNames = glob.glob("*." + fileExtension)
        for fileName in fileNames:
            os.remove(fileName)


def runExample(exampleName, mocker):
    try:
        mockPlotShow = mocker.patch("matplotlib.pyplot.show")
        importlib.import_module(exampleName)
    finally:
        cleanExampleFiles()


def runExampleWithError(exampleName, mocker, expectedErrorMessage):
    try:
        mockPlotShow = mocker.patch("matplotlib.pyplot.show")
        with pytest.raises(SystemExit) as excinfo:
            importlib.import_module(exampleName)
        assert excinfo.value.code == expectedErrorMessage
    finally:
        cleanExampleFiles()


def test_fixedAdrQstls(mocker):
    runExampleWithError(
        "fixedAdrQstls", mocker, "Error while solving the dielectric theory"
    )


def test_fixedAdrQstlsIet(mocker):
    runExample("fixedAdrQstlsIet", mocker)


def test_fixedAdrQVSStls(mocker):
    runExample("fixedAdrQVSStls", mocker)


def test_initialGuessQstls(mocker):
    runExample("initialGuessQstls", mocker)


def test_initialGuessQstlsIet(mocker):
    runExample("initialGuessQstlsIet", mocker)


def test_initialGuessStls(mocker):
    runExample("initialGuessStls", mocker)


def test_solveQuantumSchemes(mocker):
    runExample("solveQuantumSchemes", mocker)


def test_solveQVSStls(mocker):
    runExample("solveQVSStls", mocker)


def test_solveRpaAndESA(mocker):
    runExample("solveRpaAndESA", mocker)


def test_solveStls(mocker):
    runExample("solveStls", mocker)


def test_solveStlsIet(mocker):
    runExample("solveStlsIet", mocker)


def test_solveVSStls(mocker):
    runExample("solveVSStls", mocker)
