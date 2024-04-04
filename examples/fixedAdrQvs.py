import qupled.quantum as qpq

# Define a QVSStls object to solve the QVSSTLS scheme
qvsstls = qpq.QVSStls(0.1,
                    1.0,
                    mixing = 0.5,
                    resolution = 0.1,
                    cutoff = 10,
                    matsubara = 200,
                    threads = 16)

# Solve the QVSSTLS scheme
qvsstls.compute()
