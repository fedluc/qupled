import os
import pytest
from qupled.classic import ESA


@pytest.fixture
def esa():
    return ESA()


@pytest.fixture
def esa_input():
    return ESA.Input(1.0, 1.0)


def test_default(esa):
    assert esa.hdfFileName == None


def test_compute(esa, esa_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.ESA.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.ESA._checkStatusAndClean")
    mockSave = mocker.patch("qupled.classic.ESA._save")
    esa.compute(esa_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1
