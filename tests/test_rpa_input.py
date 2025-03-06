import math
import pytest
from qupled.native import RpaInput


@pytest.fixture
def rpa_input_instance():
    return RpaInput()


def test_init(rpa_input_instance):
    assert hasattr(rpa_input_instance, "coupling")
    assert hasattr(rpa_input_instance, "degeneracy")
    assert hasattr(rpa_input_instance, "theory")
    assert hasattr(rpa_input_instance, "integral_error")
    assert hasattr(rpa_input_instance, "integral_strategy")
    assert hasattr(rpa_input_instance, "threads")
    assert hasattr(rpa_input_instance, "chemical_potential")
    assert hasattr(rpa_input_instance, "matsubara")
    assert hasattr(rpa_input_instance, "resolution")
    assert hasattr(rpa_input_instance, "cutoff")
    assert hasattr(rpa_input_instance, "frequency_cutoff")


def test_defaults(rpa_input_instance):
    assert math.isnan(rpa_input_instance.coupling)
    assert math.isnan(rpa_input_instance.degeneracy)
    assert rpa_input_instance.theory == ""
    assert math.isnan(rpa_input_instance.integral_error)
    assert rpa_input_instance.integral_strategy == ""
    assert rpa_input_instance.threads == 0
    assert rpa_input_instance.chemical_potential.size == 0
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
    allowed_theories = [
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
    for theory in allowed_theories:
        rpa_input_instance.theory = theory
        this_theory = rpa_input_instance.theory
        assert this_theory == theory
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.theory = "dummyTheory"
    assert excinfo.value.args[0] == "Invalid dielectric theory: dummyTheory"


def test_integral_error(rpa_input_instance):
    rpa_input_instance.integral_error = 1.0
    integral_error = rpa_input_instance.integral_error
    assert integral_error == 1.0
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.integral_error = 0.0
    assert (
        excinfo.value.args[0]
        == "The accuracy for the integral computations must be larger than zero"
    )


def test_integral_strategy(rpa_input_instance):
    allowed_strategies = ["full", "segregated"]
    for strategy in allowed_strategies:
        rpa_input_instance.integral_strategy = strategy
        this_scheme = rpa_input_instance.integral_strategy
        assert this_scheme == strategy
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.integral_strategy = "dummyScheme"
    assert excinfo.value.args[0] == "Unknown scheme for 2D integrals: dummyScheme"


def test_threads(rpa_input_instance):
    rpa_input_instance.threads = 1
    threads = rpa_input_instance.threads
    assert threads == 1
    with pytest.raises(RuntimeError) as excinfo:
        rpa_input_instance.threads = 0
    assert excinfo.value.args[0] == "The number of threads must be larger than zero"


def test_chemical_potential(rpa_input_instance):
    rpa_input_instance.chemical_potential = [-10, 10]
    chemical_potential = rpa_input_instance.chemical_potential
    assert all(x == y for x, y in zip(chemical_potential, [-10, 10]))
    for cp in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            rpa_input_instance.chemical_potential = cp
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


def test_is_equal_default(rpa_input_instance):
    assert not rpa_input_instance.is_equal(rpa_input_instance)


def test_is_equal(rpa_input_instance):
    this_rpa = RpaInput()
    this_rpa.coupling = 2.0
    this_rpa.degeneracy = 1.0
    this_rpa.integral_error = 0.1
    this_rpa.threads = 1
    this_rpa.theory = "STLS"
    this_rpa.matsubara = 1
    this_rpa.resolution = 0.1
    this_rpa.cutoff = 1.0
    this_rpa.error = 0.1
    assert this_rpa.is_equal(this_rpa)


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
    assert "Frequency cutoff = nan" in captured
