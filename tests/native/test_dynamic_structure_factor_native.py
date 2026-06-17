import numpy as np
import pytest

from qupled.native import (
    Dimension,
    Input,
    Rpa,
    compute_dsf,
    compute_first_moment_from_dsf_adaptive,
    compute_ideal_dynamic_response,
    compute_itcf,
    compute_itcf_from_dsf_adaptive,
    compute_third_moment_from_dsf_adaptive,
)


@pytest.mark.native
def test_compute_ideal_dynamic_response():
    real, imaginary = compute_ideal_dynamic_response(
        np.array([1.2]), np.array([0.0, 0.7]), 0.8, 0.3, 1.0e-8
    )

    assert real.shape == (1, 2)
    assert imaginary.shape == (1, 2)
    assert real[0, 0] == pytest.approx(0.4803089457348037)
    assert real[0, 1] == pytest.approx(0.4245138149602567)
    assert imaginary[0, 0] == 0.0
    assert imaginary[0, 1] == pytest.approx(0.20032386853123246)


@pytest.mark.native
def test_compute_dsf_and_itcf_transform():
    inputs = Input()
    inputs.theory = "STLS"
    inputs.dimension = Dimension.D3
    inputs.coupling = 2.0
    inputs.degeneracy = 0.8
    inputs.integral_error = 1.0e-8

    frequency = np.array([0.0, 0.7, 1.4])
    dsf = compute_dsf(inputs, np.array([1.2]), frequency, 0.3, np.array([[0.25]]))
    itcf = compute_itcf_from_dsf_adaptive(
        inputs,
        np.array([1.2]),
        frequency,
        np.array([0.0, 1.25]),
        0.3,
        np.array([[0.25]]),
        dsf,
    )

    assert dsf.shape == (1, 3)
    assert dsf[0, 0] == pytest.approx(0.04345726902387644)
    assert dsf[0, 1] == pytest.approx(0.06462738121282116)
    assert itcf.shape == (1, 2)
    assert np.all(np.isfinite(itcf))
    assert np.all(itcf > 0.0)


@pytest.mark.native
def test_dsf_itcf_agrees_with_matsubara_itcf():
    inputs = Input()
    inputs.theory = "RPA"
    inputs.dimension = Dimension.D3
    inputs.coupling = 2.0
    inputs.degeneracy = 1.0
    inputs.integral_error = 1.0e-6
    inputs.integral_strategy = "full"
    inputs.chemical_potential = [-10.0, 10.0]
    inputs.cutoff = 10.0
    inputs.resolution = 1.0
    inputs.frequency_cutoff = 40.0
    inputs.matsubara = 512
    inputs.threads = 1
    scheme = Rpa(inputs)
    assert scheme.compute() == 0

    indexes = [0, 1, 2, 3, 10]
    wvg = scheme.wvg[indexes]
    lfc = scheme.lfc[indexes]
    idr = scheme.idr[indexes]
    frequency = np.arange(0.0, 40.01, 0.05)
    tau = np.array([0.5])
    dsf = compute_dsf(inputs, wvg, frequency, scheme.chemical_potential, lfc)
    spectral_itcf = compute_itcf_from_dsf_adaptive(
        inputs, wvg, frequency, tau, scheme.chemical_potential, lfc, dsf
    )
    matsubara_itcf = compute_itcf(
        inputs,
        wvg,
        tau,
        scheme.chemical_potential,
        idr,
        lfc,
    )

    assert spectral_itcf == pytest.approx(matsubara_itcf, abs=1.0e-5)


@pytest.mark.native
def test_adaptive_dsf_transform_resolves_narrow_rpa_peak():
    inputs = Input()
    inputs.theory = "RPA"
    inputs.dimension = Dimension.D3
    inputs.coupling = 2.0
    inputs.degeneracy = 1.0
    inputs.integral_error = 1.0e-6
    inputs.integral_strategy = "full"
    inputs.chemical_potential = [-10.0, 10.0]
    inputs.cutoff = 10.0
    inputs.resolution = 0.01
    inputs.frequency_cutoff = 3.0
    inputs.matsubara = 512
    inputs.threads = 1
    scheme = Rpa(inputs)
    assert scheme.compute() == 0

    indexes = [24]
    frequency = np.arange(0.0, 3.01, 0.01)
    dsf = compute_dsf(
        inputs,
        scheme.wvg[indexes],
        frequency,
        scheme.chemical_potential,
        scheme.lfc[indexes],
    )
    adaptive = compute_itcf_from_dsf_adaptive(
        inputs,
        scheme.wvg[indexes],
        frequency,
        np.array([0.0]),
        scheme.chemical_potential,
        scheme.lfc[indexes],
        dsf,
    )
    third_moment = compute_third_moment_from_dsf_adaptive(
        inputs,
        scheme.wvg[indexes],
        frequency,
        scheme.chemical_potential,
        scheme.lfc[indexes],
        dsf,
    )
    coarse_frequency = np.arange(0.0, 3.01, 0.05)
    coarse_dsf = compute_dsf(
        inputs,
        scheme.wvg[indexes],
        coarse_frequency,
        scheme.chemical_potential,
        scheme.lfc[indexes],
    )
    coarse_third_moment = compute_third_moment_from_dsf_adaptive(
        inputs,
        scheme.wvg[indexes],
        coarse_frequency,
        scheme.chemical_potential,
        scheme.lfc[indexes],
        coarse_dsf,
    )

    assert adaptive[:, 0] == pytest.approx(scheme.ssf[indexes], abs=1.0e-3)
    assert third_moment.shape == (len(indexes),)
    assert np.all(np.isfinite(third_moment))
    assert np.all(third_moment > 0.0)
    assert coarse_third_moment == pytest.approx(third_moment, rel=8.0e-2)


@pytest.mark.native
@pytest.mark.parametrize("theory", ["RPA", "STLS"])
def test_adaptive_first_moment_satisfies_f_sum_rule(theory):
    inputs = Input()
    inputs.theory = theory
    inputs.dimension = Dimension.D3
    inputs.coupling = 3.23
    inputs.degeneracy = 1.0
    inputs.integral_error = 1.0e-6
    inputs.integral_strategy = "full"
    inputs.chemical_potential = [-10.0, 10.0]
    inputs.cutoff = 4.0
    inputs.resolution = 0.1
    inputs.frequency_cutoff = 50.0
    inputs.matsubara = 256
    inputs.threads = 1
    scheme = Rpa(inputs)
    assert scheme.compute() == 0

    indexes = [5, 10, 20, 30]
    wvg = scheme.wvg[indexes]
    lfc = scheme.lfc[indexes]
    if theory == "STLS":
        lfc = 0.25 * np.ones_like(lfc)
    frequency = np.arange(0.0, 50.01, 0.05)
    dsf = compute_dsf(
        inputs, wvg, frequency, scheme.chemical_potential, lfc
    )
    first_moment = compute_first_moment_from_dsf_adaptive(
        inputs, wvg, frequency, scheme.chemical_potential, lfc, dsf
    )

    assert first_moment == pytest.approx(wvg**2, rel=2.0e-3, abs=2.0e-4)
