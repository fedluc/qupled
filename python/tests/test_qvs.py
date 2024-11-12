import os
import pytest
import numpy as np
from qupled.qupled import QVSStls as QVSStlsNative
from qupled.util import Hdf
from qupled.quantum import QVSStls


@pytest.fixture
def qvsstls():
    return QVSStls()


@pytest.fixture
def qvsstls_input():
    return QVSStls.Input(1.0, 1.0)


def test_default(qvsstls):
    assert qvsstls.hdfFileName is None


def test_compute(qvsstls, qvsstls_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockUnpack = mocker.patch("qupled.quantum.QVSStls._unpackFixedAdrFiles")
    mockCompute = mocker.patch("qupled.quantum.QVSStls._compute")
    mockSave = mocker.patch("qupled.quantum.QVSStls._save")
    mockZip = mocker.patch("qupled.quantum.QVSStls._zipFixedAdrFiles")
    mockClean = mocker.patch("qupled.quantum.QVSStls._cleanFixedAdrFiles")
    qvsstls.compute(qvsstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockUnpack.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1
    assert mockZip.call_count == 1
    assert mockClean.call_count == 1


def test_unpackFixedAdrFiles_no_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockExtractAll = mocker.patch("qupled.quantum.zf.ZipFile.extractall")
    qvsstls._unpackFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 0
    assert mockExtractAll.call_count == 0


def test_unpackFixedAdrFiles_with_files(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockExtractAll = mocker.patch("qupled.quantum.zf.ZipFile.extractall")
    qvsstls_input.fixed = "testFile.zip"
    qvsstls._unpackFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockExtractAll.call_count == 1


def test_zipFixedAdrFiles_no_file(qvsstls, qvsstls_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockGlob = mocker.patch(
        "qupled.quantum.glob", return_value={"binFile1", "binFile2"}
    )
    mockRemove = mocker.patch("os.remove")
    mockWrite = mocker.patch("qupled.quantum.zf.ZipFile.write")
    qvsstls._zipFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockGlob.call_count == 1
    assert mockRemove.call_count == 2
    assert mockWrite.call_count == 2


def test_cleanFixedAdrFiles_no_files(qvsstls, qvsstls_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockRemove = mocker.patch("qupled.quantum.rmtree")
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 0


def test_cleanFixedAdrFiles_with_files(qvsstls, qvsstls_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockIsDir = mocker.patch("os.path.isdir", return_value=True)
    mockRemove = mocker.patch("qupled.quantum.rmtree")
    qvsstls._cleanFixedAdrFiles(qvsstls_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 1


def test_save(qvsstls, qvsstls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = QVSStlsNative(qvsstls_input.toNative())
        qvsstls.hdfFileName = qvsstls._getHdfFile(scheme.inputs)
        qvsstls._save(scheme)
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qvsstls.hdfFileName)
        inspectData = Hdf().inspect(qvsstls.hdfFileName)
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
            "adr",
            "alpha",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qvsstls.hdfFileName)


def test_setFreeEnergyIntegrand(mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"fxcGrid": arr1D, "fxci": arr2D, "alpha": arr1D},
    )
    fxc = QVSStls.getFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(fxc.grid, arr1D)
    assert np.array_equal(fxc.alpha, arr1D)
    assert np.array_equal(fxc.integrand, arr2D)
