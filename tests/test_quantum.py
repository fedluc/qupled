from qupled.quantum import Qstls, QstlsIet, QVSStls


def test_qstls_initialization():
    qstls_instance = Qstls()
    assert isinstance(qstls_instance, Qstls)


def test_qstlsiet_initialization():
    qstlsiet_instance = QstlsIet()
    assert isinstance(qstlsiet_instance, QstlsIet)


def test_qvsstls_initialization():
    qvsstls_instance = QVSStls()
    assert isinstance(qvsstls_instance, QVSStls)
