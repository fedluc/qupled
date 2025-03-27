import os
import glob
from qupled.native import Qstls, Stls, QstlsInput


def test_qstls_properties():
    assert issubclass(Qstls, Stls)
    scheme = Qstls(QstlsInput())
    assert hasattr(scheme, "adr")


def test_qstls_compute():
    inputs = QstlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "QSTLS"
    inputs.chemical_potential = [-10, 10]
    inputs.cutoff = 5.0
    inputs.matsubara = 32
    inputs.resolution = 0.1
    inputs.integral_error = 1.0e-5
    inputs.threads = 16
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.output_frequency = 2
    scheme = Qstls(inputs)
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
        assert scheme.recovery == "recovery_rs1.000_theta1.000_QSTLS.bin"
        assert os.path.isfile(scheme.recovery)
        assert scheme.rdf(scheme.wvg).size == nx
    finally:
        fixed_file = "adr_fixed_theta1.000_matsubara32.bin"
        if os.path.isfile(scheme.recovery):
            os.remove(scheme.recovery)
        if os.path.isfile(fixed_file):
            os.remove(fixed_file)


def test_qstls_iet_compute():
    iet_schemes = {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}
    for scheme_name in iet_schemes:
        inputs = QstlsInput()
        inputs.coupling = 10.0
        inputs.degeneracy = 1.0
        inputs.theory = scheme_name
        inputs.chemical_potential = [-10, 10]
        inputs.cutoff = 5.0
        inputs.matsubara = 16
        inputs.resolution = 0.1
        inputs.integral_error = 1.0e-5
        inputs.threads = 16
        inputs.error = 1.0e-5
        inputs.mixing = 0.5
        inputs.iterations = 1000
        inputs.output_frequency = 2
        scheme = Qstls(inputs)
        scheme.compute()
        try:
            nx = scheme.wvg.size
            assert nx >= 3
            assert scheme.idr.shape[0] == nx
            assert scheme.idr.shape[1] == inputs.matsubara
            assert scheme.sdr.size == nx
            assert scheme.slfc.size == nx
            assert scheme.ssf.size == nx
            assert scheme.ssf_HF.size == nx
            recovery = "recovery_rs10.000_theta1.000_" + scheme_name + ".bin"
            assert scheme.recovery == recovery
            assert os.path.isfile(scheme.recovery)
            assert scheme.rdf(scheme.wvg).size == nx
        finally:
            if os.path.isfile(scheme.recovery):
                os.remove(scheme.recovery)
            file_names = glob.glob("adr_fixed*.bin")
            for file_name in file_names:
                os.remove(file_name)
