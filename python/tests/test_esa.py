import os
import pytest
import set_path
from qupled.classic import ESA
from qupled.classic import ClassicScheme


@pytest.fixture
def esa_instance():
    return ESA(1.0, 1.0)


def test_default(esa_instance):
    assert issubclass(ESA, ClassicScheme)
    assert all(x == y for x, y in zip(esa_instance.allowedTheories, ["ESA"]))
    assert esa_instance.inputs.coupling == 1.0
    assert esa_instance.inputs.degeneracy == 1.0
    assert esa_instance.inputs.theory == "ESA"
    assert all(
        x == y for x, y in zip(esa_instance.inputs.chemicalPotential, [-10.0, 10.0])
    )
    assert esa_instance.inputs.cutoff == 10.0
    assert esa_instance.inputs.matsubara == 128
    assert esa_instance.inputs.resolution == 0.1
    assert esa_instance.inputs.intError == 1.0e-5
    assert esa_instance.inputs.threads == 1
    assert esa_instance.scheme is None
    assert esa_instance.hdfFileName is None


def test_set_input():
    esa = ESA(2.0, 0.5, [-5, 5], 20, 32, 0.01)
    assert esa.inputs.coupling == 2.0
    assert esa.inputs.degeneracy == 0.5
    assert esa.inputs.theory == "ESA"
    assert all(x == y for x, y in zip(esa.inputs.chemicalPotential, [-5.0, 5.0]))
    assert esa.inputs.cutoff == 20
    assert esa.inputs.matsubara == 32
    assert esa.inputs.resolution == 0.01


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
