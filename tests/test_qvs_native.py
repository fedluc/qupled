import os
from qupled.native import QVSStls, QVSStlsInput, Rpa


def test_qvsstls_properties():
    assert issubclass(QVSStls, Rpa)
    inputs = QVSStlsInput()
    inputs.coupling = 1.0
    inputs.coupling_resolution = 0.1
    scheme = QVSStls(inputs)
    assert hasattr(scheme, "free_energy_integrand")
    assert hasattr(scheme, "free_energy_grid")
    assert hasattr(scheme, "adr")


def test_qvsstls_compute():
    inputs = QVSStlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "QVSSTLS"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 5.0
    inputs.matsubara = 32
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 16
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.output_frequency = 10
    inputs.coupling_resolution = 0.1
    inputs.degeneracy_resolution = 0.1
    inputs.error_alpha = 1.0e-3
    inputs.iterations_alpha = 50
    inputs.alpha = [0.5, 1.0]
    scheme = QVSStls(inputs)
    scheme.compute()
    try:
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.adr.shape[0] == nx
        assert scheme.adr.shape[1] == inputs.matsubara
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.ssf_HF.size == nx
        assert scheme.recovery == "recovery_rs1.000_theta1.000_QVSSTLS.bin"
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        fixed_filem = "THETA_DOWN.bin"
        fixed_file = "THETA.bin"
        fixed_filep = "THETA_UP.bin"
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recovery)
        if os.path.isfile(fixed_filem):
            os.remove(fixed_filem)
        if os.path.isfile(fixed_file):
            os.remove(fixed_file)
        if os.path.isfile(fixed_filep):
            os.remove(fixed_filep)
