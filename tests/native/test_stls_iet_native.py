from qupled.native import Stls, StlsIet, StlsIetInput


def test_stls_properties():
    assert issubclass(StlsIet, Stls)
    scheme = StlsIet(StlsIetInput())
    assert hasattr(scheme, "bf")


def test_stls_iet_compute():
    iet_schemes = {"STLS-HNC", "STLS-IOI", "STLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = StlsIetInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = scheme_name
        inputs.chemical_potential = [-10, 10]
        inputs.cutoff = 5.0
        inputs.matsubara = 32
        inputs.resolution = 0.1
        inputs.integral_error = 1.0e-5
        inputs.threads = 1
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        scheme = StlsIet(inputs)
        scheme.compute()
        nx = scheme.wvg.size
        assert nx >= 3
        assert scheme.idr.shape[0] == nx
        assert scheme.idr.shape[1] == inputs.matsubara
        assert scheme.sdr.size == nx
        assert scheme.slfc.size == nx
        assert scheme.ssf.size == nx
        assert scheme.bf.size == nx
        assert scheme.rdf(scheme.wvg).size == nx
