import os
import pytest
import math
import set_path
import qupled.qupled as qp
import qupled.classic as qpc

def tolerance():
    return 1e-10

def test_vsstls_properties():
    issubclass(qp.VSStls, qp.Rpa)
    inputs = qpc.VSStls(1.0, 1.0).inputs
    scheme = qp.VSStls(inputs)
    assert hasattr(scheme, "freeEnergyIntegrand")
    assert hasattr(scheme, "freeEnergyGrid")
    
def test_vsstls_compute():
    inputs = qpc.VSStls(1.0, 1.0,
                        couplingResolution=0.1,
                        degeneracyResolution=0.1,
                        cutoff=5).inputs
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
        assert scheme.recovery == ""
        assert scheme.rdf(scheme.wvg).size == nx
        assert math.isclose(scheme.uInt, -0.5223739856, rel_tol=tolerance())
    finally:
        if (os.path.isfile(scheme.recovery)) : os.remove(scheme.recover)

