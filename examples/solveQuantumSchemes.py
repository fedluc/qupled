import numpy as np
import qupled.quantum as qpq

# Define a Qstls object to solve the QSTLS scheme
qstls = qpq.Qstls(10.0, 1.0,
                  mixing = 0.5,
                  resolution = 0.1,
                  cutoff = 10,
                  matsubara = 16,
                  threads = 16)

# Solve the QSTLS scheme
qstls.compute()

# Plot the density responses and the static local field correction 
qstls.plot(["idr", "adr"], matsubara = np.arange(1, 10, 2))

# Define a QstlsIet object to solve one of the QSTLS-IET schemes
qstls = qpq.QstlsIet(30.0, 1.0,
                     "QSTLS-LCT",
                     mixing = 0.2,
                     resolution = 0.1,
                     cutoff = 10,
                     matsubara = 16,
                     scheme2DIntegrals = "segregated",
                     threads = 16)

# solve the QSTLS-IET scheme
qstls.compute()
