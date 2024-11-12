import os
import pytest
import qupled.qupled as qp


def test_vsstls_properties():
    assert issubclass(qp.VSStls, qp.Rpa)
    inputs = qp.VSStlsInput()
    inputs.coupling = 1.0
    inputs.couplingResolution = 0.1
    scheme = qp.VSStls(inputs)
    assert hasattr(scheme, "freeEnergyIntegrand")
    assert hasattr(scheme, "freeEnergyGrid")


def test_vsstls_compute():
    inputs = qp.VSStlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "VSSTLS"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 5.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.outputFrequency = 10
    inputs.couplingResolution = 0.1
    inputs.degeneracyResolution = 0.1
    inputs.errorAlpha = 1.0e-3
    inputs.iterationsAlpha = 50
    inputs.alpha = [0.5, 1.0]
    scheme = qp.VSStls(inputs)
    scheme.compute()
    try:
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.ssfHF.size == nx
        assert scheme.recovery == "recovery_rs1.000_theta1.000_VSSTLS.bin"
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recover)
