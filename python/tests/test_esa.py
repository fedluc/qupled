import os
import pytest
import set_path
import qupled.qupled as qp
from qupled.classic import ESA, ClassicSchemeNew


@pytest.fixture
def esa_instance():
    assert issubclass(qp.Stls, qp.Rpa)
    inputs = qp.RpaInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "ESA"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    return ESA(inputs)


def test_default(esa_instance):
    assert issubclass(ESA, ClassicSchemeNew)
    assert issubclass(ESA, qp.ESA)
    assert all(x == y for x, y in zip(esa_instance.allowedTheories, ["ESA"]))
    assert esa_instance.hdfFileName is None


def test_compute(esa_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.ESA._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.ESA.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.ESA._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.ESA._setHdfFile")
    mockSave = mocker.patch("qupled.classic.ESA._save")
    esa_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1
