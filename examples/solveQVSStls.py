import numpy as np
import qupled.quantum as qpq

# Define a Qstls object to solve the QSTLS scheme
qvsstls = qpq.QVSStls(10.0, 1.0,
                  mixing = 0.5,
                  resolution = 0.1,
                  cutoff = 10,
                  matsubara = 16,
                  threads = 1)

# Solve the QSTLS scheme
qvsstls.compute()

# Plot the density responses and the static local field correction 
qvsstls.plot(["idr", "adr", "slfc"], matsubara = np.arange(1, 10, 2))

