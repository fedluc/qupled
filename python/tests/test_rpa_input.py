import os
import math
import pytest
import qupled.qupled as qp


@pytest.fixture
def rpa_input_instance():
    return qp.RpaInput()


def test_init(rpa_input_instance):
    assert hasattr(rpa_input_instance, "coupling")
    assert hasattr(rpa_input_instance, "degeneracy")
    assert hasattr(rpa_input_instance, "theory")
    assert hasattr(rpa_input_instance, "intError")
    assert hasattr(rpa_input_instance, "int2DScheme")
    assert hasattr(rpa_input_instance, "threads")
    assert hasattr(rpa_input_instance, "chemicalPotential")
    assert hasattr(rpa_input_instance, "matsubara")
    assert hasattr(rpa_input_instance, "resolution")
    assert hasattr(rpa_input_instance, "cutoff")


def test_defaults(rpa_input_instance):
    assert math.isnan(rpa_input_instance.coupling)
    assert math.isnan(rpa_input_instance.degeneracy)
    assert rpa_input_instance.theory == ""
    assert math.isnan(rpa_input_instance.intError)
    assert rpa_input_instance.int2DScheme == ""
    assert rpa_input_instance.threads == 0
    assert rpa_input_instance.chemicalPotential.size == 0
    assert rpa_input_instance.matsubara == 0
    assert math.isnan(rpa_input_instance.resolution)
    assert math.isnan(rpa_input_instance.cutoff)


def test_coupling(rpa_input_instance):
    rpa_input_instance.coupling = 1.0
    coupling = rpa_input_instance.coupling
    assert coupling == 1.0
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.coupling = -1.0
    assert excinfo.value.args[0] == "The quantum coupling parameter can't be negative"


def test_degeneracy(rpa_input_instance):
    rpa_input_instance.degeneracy = 1.0
    degeneracy = rpa_input_instance.degeneracy
    assert degeneracy == 1.0
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.degeneracy = -1.0
    assert excinfo.value.args[0] == "The quantum degeneracy parameter can't be negative"


def test_theory(rpa_input_instance):
    allowedTheories = [
        "RPA",
        "ESA",
        "STLS",
        "STLS-HNC",
        "STLS-IOI",
        "STLS-LCT",
        "VSSTLS",
        "QSTLS",
        "QSTLS-HNC",
        "QSTLS-IOI",
        "QSTLS-LCT",
    ]
    for theory in allowedTheories:
        rpa_input_instance.theory = theory
        thisTheory = rpa_input_instance.theory
        assert thisTheory == theory
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.theory = "dummyTheory"
    assert excinfo.value.args[0] == "Invalid dielectric theory: dummyTheory"


def test_intError(rpa_input_instance):
    rpa_input_instance.intError = 1.0
    intError = rpa_input_instance.intError
    assert intError == 1.0
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.intError = 0.0
    assert (
        excinfo.value.args[0]
        == "The accuracy for the integral computations must be larger than zero"
    )


def test_int2DScheme(rpa_input_instance):
    allowedSchemes = ["full", "segregated"]
    for scheme in allowedSchemes:
        rpa_input_instance.int2DScheme = scheme
        thisScheme = rpa_input_instance.int2DScheme
        assert thisScheme == scheme
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.int2DScheme = "dummyScheme"
    assert excinfo.value.args[0] == "Unknown scheme for 2D integrals: dummyScheme"


def test_threads(rpa_input_instance):
    rpa_input_instance.threads = 1
    threads = rpa_input_instance.threads
    assert threads == 1
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.threads = 0
    assert excinfo.value.args[0] == "The number of threads must be larger than zero"


def test_chemicalPotential(rpa_input_instance):
    rpa_input_instance.chemicalPotential = [-10, 10]
    chemicalPotential = rpa_input_instance.chemicalPotential
    assert all(x == y for x, y in zip(chemicalPotential, [-10, 10]))
    for cp in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            rpa_input_instance.chemicalPotential = cp
        assert (
            excinfo.value.args[0] == "Invalid guess for chemical potential calculation"
        )


def test_matsubara(rpa_input_instance):
    rpa_input_instance.matsubara = 1
    matsubara = rpa_input_instance.matsubara
    assert matsubara == 1
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.matsubara = -1
    assert (
        excinfo.value.args[0] == "The number of matsubara frequencies can't be negative"
    )


def test_resolution(rpa_input_instance):
    rpa_input_instance.resolution = 0.01
    resolution = rpa_input_instance.resolution
    assert resolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.resolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The wave-vector grid resolution must be larger than zero"
    )


def test_cutoff(rpa_input_instance):
    rpa_input_instance.cutoff = 0.01
    cutoff = rpa_input_instance.cutoff
    assert cutoff == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.cutoff = -0.1
    assert (
        excinfo.value.args[0] == "The wave-vector grid cutoff must be larger than zero"
    )


def test_isEqual_default(rpa_input_instance):
    assert not rpa_input_instance.isEqual(rpa_input_instance)


def test_isEqual(rpa_input_instance):
    thisRpa = qp.RpaInput()
    thisRpa.coupling = 2.0
    thisRpa.degeneracy = 1.0
    thisRpa.intError = 0.1
    thisRpa.threads = 1
    thisRpa.theory = "STLS"
    thisRpa.matsubara = 1
    thisRpa.resolution = 0.1
    thisRpa.cutoff = 1.0
    thisRpa.error = 0.1
    assert thisRpa.isEqual(thisRpa)


def test_print(rpa_input_instance, capfd):
    rpa_input_instance.print()
    captured = capfd.readouterr().out
    captured = captured.split("\n")
    assert "Coupling parameter = nan" in captured
    assert "Degeneracy parameter = nan" in captured
    assert "Number of OMP threads = 0" in captured
    assert "Scheme for 2D integrals = " in captured
    assert "Integral relative error = nan" in captured
    assert "Theory to be solved = " in captured
    assert "Guess for chemical potential = " in captured
    assert "Number of Matsubara frequencies = 0" in captured
    assert "Wave-vector resolution = nan" in captured
    assert "Wave-vector cutoff = nan" in captured
