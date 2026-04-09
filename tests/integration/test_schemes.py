import os
import shutil

import pytest

from qupled.database.database_handler import DataBaseHandler, DATABASE_DIRECTORY


@pytest.fixture(autouse=True)
def run_after_each_test():
    yield
    output_dir = DATABASE_DIRECTORY
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)


def assert_internal_energy(expected_internal_energy, tolerance=1e-5):
    db_handler = DataBaseHandler()
    for run_id, expected_uint in expected_internal_energy.items():
        uint = db_handler.scheme_tables.get_results(run_id)["uint"]
        assert abs(uint - expected_uint) < tolerance * abs(expected_uint)


@pytest.mark.integration
def test_fixed_adr_qstls():
    expected_internal_energy = {
        1: -0.06915828119355959,
        2: -0.06915828119355959,
        3: -0.036272826904439975,
        4: -0.03545665197810883,
    }
    from qupled.schemes import qstls

    # Solve qstls
    scheme = qstls.Solver()
    inputs = qstls.Input(10.0, 1.0, mixing=0.5, matsubara=16, threads=16)
    scheme.compute(inputs)
    # Solve qslts with the same inputs (this will reuse the fixed component from the previous computation)
    scheme.compute(inputs)
    # Change the coupling parameter
    inputs.coupling = 20.0
    # Compute with the updated coupling parameter (this will recompute the fixed component)
    scheme.compute(inputs)
    # Change the degeneracy parameter
    inputs.degeneracy = 2.0
    # Compute with the update degeneracy parameter (this will recompute the fixed component)
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_initial_guess_stls():
    expected_internal_energy = {1: -0.0696211823813233, 2: -0.06962120294563244}
    from qupled.schemes import stls

    # Define the object used to solve the scheme
    scheme = stls.Solver()
    # Define the input parameters
    inputs = stls.Input(10.0, 1.0, mixing=0.2)
    # Solve scheme
    scheme.compute(inputs)
    # Create a custom initial guess from the output files of the previous run
    inputs.guess = scheme.get_initial_guess(scheme.run_id)
    # Solve the scheme again with the new initial guess
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_esa():
    expected_internal_energy = {1: -0.0700592697600796}
    from qupled.schemes import esa

    scheme = esa.Solver()
    scheme.compute(esa.Input(10.0, 1.0))
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_rpa():
    expected_internal_energy = {1: -0.09410276327343531}
    from qupled.schemes import rpa

    scheme = rpa.Solver()
    scheme.compute(rpa.Input(10.0, 1.0))
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_qstls():
    expected_internal_energy = {1: -0.06915828119355959}
    from qupled.schemes import qstls

    # Define a Qstls object to solve the QSTLS scheme
    scheme = qstls.Solver()
    # Define the input parameters
    inputs = qstls.Input(10.0, 1.0, mixing=0.5, matsubara=16, threads=16)
    # Solve the QSTLS scheme
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_qstls_iet():
    expected_internal_energy = {1: -0.07152825527327064}
    from qupled.schemes import qstlsiet

    # Define a QstlsIet object to solve the QSTLS-IET scheme
    scheme = qstlsiet.Solver()
    # Define the input parameters for one of the QSTLS-IET schemes
    inputs = qstlsiet.Input(
        10.0,
        1.0,
        theory="QSTLS-LCT",
        mixing=0.5,
        matsubara=16,
        threads=16,
        integral_strategy="segregated",
    )
    # solve the QSTLS-IET scheme
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_qvsstls():
    expected_internal_energy = {8: -0.282671921632354}
    from qupled.schemes import qvsstls

    # Define the object used to solve the scheme
    scheme = qvsstls.Solver()
    # Define the input parameters
    inputs = qvsstls.Input(
        1.0,
        1.0,
        mixing=0.5,
        matsubara=16,
        alpha=[-0.2, 0.4],
        iterations=100,
        threads=16,
    )
    # Solve scheme for rs = 1.0
    scheme.compute(inputs)
    # Load the free energy integrand computed for rs = 1.0
    fxci = scheme.get_free_energy_integrand(scheme.run_id)
    # Setup a new  simulation for rs=2.0
    inputs.coupling = 2.0
    inputs.alpha = [0.1, 0.5]
    inputs.free_energy_integrand = fxci
    # Solve scheme for rs = 2.0
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_stls():
    expected_internal_energy = {1: -0.0696212348180385}
    from qupled.schemes import stls

    scheme = stls.Solver()
    inputs = stls.Input(10.0, 1.0, mixing=0.5)
    scheme.compute(inputs)
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_vsstls():
    expected_internal_energy = {19: -0.1333040908470305}
    from qupled.schemes import vsstls

    scheme = vsstls.Solver()
    inputs = vsstls.Input(2.0, 1.0, mixing=0.5, alpha=[-0.2, 0.2])
    scheme.compute(inputs)
    # Load the free energy integrand computed for rs = 2.0
    fxci = scheme.get_free_energy_integrand(scheme.run_id)
    # Setup a new VSStls simulation for rs = 5.0
    inputs.coupling = 5.0
    inputs.alpha = [0.5, 0.7]
    inputs.free_energy_integrand = fxci
    # Solve again
    scheme.compute(inputs)
    # Assert
    assert_internal_energy(expected_internal_energy)

@pytest.mark.integration
def test_solve_2D_ground():
    expected_internal_energy = {1: -0.060018666228726536, 2: -0.1607977129720737, 3: -0.09514934062352233}
    from qupled.schemes import rpa, stls, hf
    from qupled.util.dimension import Dimension

    solver = hf.Solver()
    solver.compute(hf.Input(10.0, 0.0, dimension=Dimension._2D))
    solver = rpa.Solver()
    solver.compute(rpa.Input(10.0, 0.0, dimension=Dimension._2D))
    solver = stls.Solver()
    solver.compute(stls.Input(10.0, 0.0, dimension=Dimension._2D, mixing=0.4))
    # Assert
    assert_internal_energy(expected_internal_energy)


@pytest.mark.integration
def test_solve_2D_finite_temperature():
    expected_internal_energy = {1: -0.03831335655062619, 2: -0.18725913760698018, 3: -0.09591039535710567}
    from qupled.schemes import rpa, stls, hf
    from qupled.util.dimension import Dimension

    solver = hf.Solver()
    solver.compute(hf.Input(10.0, 1.0, dimension=Dimension._2D))
    solver = rpa.Solver()
    solver.compute(rpa.Input(10.0, 1.0, dimension=Dimension._2D))
    solver = stls.Solver()
    solver.compute(stls.Input(10.0, 1.0, dimension=Dimension._2D, mixing=0.3))
    # Assert
    assert_internal_energy(expected_internal_energy)