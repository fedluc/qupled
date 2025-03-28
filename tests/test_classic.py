from qupled.classic import ESA, Rpa, Stls, StlsIet, VSStls


def test_esa_initialization():
    esa = ESA()
    assert isinstance(esa, ESA)


def test_rpa_initialization():
    rpa = Rpa()
    assert isinstance(rpa, Rpa)


def test_stls_initialization():
    stls = Stls()
    assert isinstance(stls, Stls)


def test_stlsiet_initialization():
    stlsiet = StlsIet()
    assert isinstance(stlsiet, StlsIet)


def test_vsstls_initialization():
    vsstls = VSStls()
    assert isinstance(vsstls, VSStls)
