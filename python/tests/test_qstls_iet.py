import os
import pytest
import numpy as np
from qupled.qupled import Qstls as QstlsNative
from qupled.util import Hdf
from qupled.quantum import QstlsIet


@pytest.fixture
def qstls_iet():
    return QstlsIet()


@pytest.fixture
def qstls_iet_input():
    return QstlsIet.Input(1.0, 1.0, "QSTLS-HNC")


def test_default(qstls_iet):
    assert qstls_iet.hdfFileName is None


def test_compute(qstls_iet, qstls_iet_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockUnpack = mocker.patch("qupled.quantum.QstlsIet._unpackFixedAdrFiles")
    mockCompute = mocker.patch("qupled.quantum.QstlsIet._compute")
    mockSave = mocker.patch("qupled.quantum.QstlsIet._save")
    mockZip = mocker.patch("qupled.quantum.QstlsIet._zipFixedAdrFiles")
    mockClean = mocker.patch("qupled.quantum.QstlsIet._cleanFixedAdrFiles")
    qstls_iet.compute(qstls_iet_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockUnpack.call_count == 1
    assert mockCompute.call_count == 1
    assert mockSave.call_count == 1
    assert mockZip.call_count == 1
    assert mockClean.call_count == 1


def test_unpackFixedAdrFiles_no_files(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockExtractAll = mocker.patch("qupled.quantum.zf.ZipFile.extractall")
    qstls_iet._unpackFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 0
    assert mockExtractAll.call_count == 0


def test_unpackFixedAdrFiles_with_files(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockExtractAll = mocker.patch("qupled.quantum.zf.ZipFile.extractall")
    qstls_iet_input.fixediet = "testFile.zip"
    qstls_iet._unpackFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockExtractAll.call_count == 1


def test_zipFixedAdrFiles_no_file(qstls_iet, qstls_iet_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockZip = mocker.patch("qupled.quantum.zf.ZipFile.__init__", return_value=None)
    mockGlob = mocker.patch(
        "qupled.quantum.glob", return_value={"binFile1", "binFile2"}
    )
    mockRemove = mocker.patch("os.remove")
    mockWrite = mocker.patch("qupled.quantum.zf.ZipFile.write")
    qstls_iet._zipFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockZip.call_count == 1
    assert mockGlob.call_count == 1
    assert mockRemove.call_count == 2
    assert mockWrite.call_count == 2


def test_cleanFixedAdrFiles_no_files(qstls_iet, qstls_iet_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockRemove = mocker.patch("qupled.quantum.rmtree")
    qstls_iet._cleanFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 0


def test_cleanFixedAdrFiles_with_files(qstls_iet, qstls_iet_input, mocker, capsys):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    mockIsDir = mocker.patch("os.path.isdir", return_value=True)
    mockRemove = mocker.patch("qupled.quantum.rmtree")
    qstls_iet._cleanFixedAdrFiles(qstls_iet_input)
    assert mockMPIIsRoot.call_count == 1
    assert mockRemove.call_count == 1


def test_save(qstls_iet, qstls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = QstlsNative(qstls_iet_input.toNative())
        qstls_iet.hdfFileName = qstls_iet._getHdfFile(scheme.inputs)
        qstls_iet._save(scheme)
        assert mockMPIIsRoot.call_count == 4
        assert os.path.isfile(qstls_iet.hdfFileName)
        inspectData = Hdf().inspect(qstls_iet.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "adr",
            "idr",
            "sdr",
            "slfc",
            "bf",
            "ssf",
            "ssfHF",
            "wvg",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qstls_iet.hdfFileName)
