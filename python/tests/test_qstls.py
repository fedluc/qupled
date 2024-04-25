import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.quantum import Qstls
from qupled.classic import Stls


@pytest.fixture
def qstls_instance():
    return Qstls(1.0, 1.0)


def test_default(qstls_instance):
    assert issubclass(Qstls, Stls)
    assert all(x == y for x, y in zip(qstls_instance.allowedTheories, ["QSTLS"]))
    assert qstls_instance.inputs.coupling == 1.0
    assert qstls_instance.inputs.degeneracy == 1.0
    assert qstls_instance.inputs.theory == "QSTLS"
    assert all(
        x == y for x, y in zip(qstls_instance.inputs.chemicalPotential, [-10.0, 10.0])
    )
    assert qstls_instance.inputs.cutoff == 10.0
    assert qstls_instance.inputs.error == 1.0e-5
    assert qstls_instance.inputs.fixed == ""
    assert qstls_instance.inputs.mixing == 1.0
    assert qstls_instance.inputs.iterations == 1000
    assert qstls_instance.inputs.matsubara == 128
    assert qstls_instance.inputs.outputFrequency == 10
    assert qstls_instance.inputs.recoveryFile == ""
    assert qstls_instance.inputs.resolution == 0.1
    assert qstls_instance.inputs.int2DScheme == "full"
    assert qstls_instance.inputs.threads == 1
    assert qstls_instance.inputs.intError == 1.0e-5
    assert qstls_instance.scheme is None
    assert qstls_instance.hdfFileName is None


def test_set_input():
    qstls = Qstls(
        2.0,
        0.5,
        [-5, 5],
        20,
        1.0e-8,
        "fixedAdrFile",
        0.5,
        None,
        100,
        32,
        100,
        "recoveryFile",
        0.01,
        "segregated",
        16,
    )
    assert qstls.inputs.coupling == 2.0
    assert qstls.inputs.degeneracy == 0.5
    assert qstls.inputs.theory == "QSTLS"
    assert all(x == y for x, y in zip(qstls.inputs.chemicalPotential, [-5.0, 5.0]))
    assert qstls.inputs.cutoff == 20
    assert qstls.inputs.error == 1.0e-8
    assert qstls.inputs.fixed == "fixedAdrFile"
    assert qstls.inputs.mixing == 0.5
    assert qstls.inputs.iterations == 100
    assert qstls.inputs.matsubara == 32
    assert qstls.inputs.outputFrequency == 100
    assert qstls.inputs.recoveryFile == "recoveryFile"
    assert qstls.inputs.resolution == 0.01
    assert qstls.inputs.int2DScheme == "segregated"
    assert qstls.inputs.threads == 16


def test_compute(qstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.quantum.Qstls._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Qstls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.quantum.Qstls._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.quantum.Qstls._setHdfFile")
    mockSave = mocker.patch("qupled.quantum.Qstls._save")
    qstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_save(qstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    qstls_instance.scheme = qp.Qstls(qstls_instance.inputs)
    qstls_instance._setHdfFile()
    try:
        qstls_instance._save()
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(qstls_instance.hdfFileName)
        inspectData = Hdf().inspect(qstls_instance.hdfFileName)
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
            "ssf",
            "ssfHF",
            "wvg",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(qstls_instance.hdfFileName)


def test_setGuess(qstls_instance, mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read", return_value={"wvg": arr, "ssf": arr}
    )
    qstls_instance.setGuess("dummyFileName")
    assert np.array_equal(qstls_instance.inputs.guess.wvg, arr)
    assert np.array_equal(qstls_instance.inputs.guess.ssf, arr)
