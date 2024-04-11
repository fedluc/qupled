import os
import pytest
import set_path
import qupled.qupled as qp
import qupled.quantum as qpq

def test_qvsstls_properties():
    assert issubclass(qp.QVSStls, qp.Rpa)
    inputs = qpq.QVSStls(1.0, 1.0).inputs
    scheme = qp.QVSStls(inputs)
    assert hasattr(scheme, "freeEnergyIntegrand")
    assert hasattr(scheme, "freeEnergyGrid")
    assert hasattr(scheme, "adr")
    
def test_qvsstls_compute():
    inputs = qpq.QVSStls(1.0, 1.0,
                         matsubara=32,
                         couplingResolution=0.1,
                         degeneracyResolution=0.1,
                         cutoff=5,
                         threads=16).inputs
    scheme = qp.QVSStls(inputs)
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
        assert scheme.ssfHF.size == nx
        assert scheme.recovery == ""
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        fixedFilem = "adr_fixed_theta0.900_matsubara32.bin"
        fixedFile = "adr_fixed_theta1.000_matsubara32.bin"
        fixedFilep = "adr_fixed_theta1.100_matsubara32.bin"
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recover)
        if (os.path.isfile(fixedFilem)) : os.remove(fixedFilem)
        if (os.path.isfile(fixedFile)) : os.remove(fixedFile)
        if (os.path.isfile(fixedFilep)) : os.remove(fixedFilep)

