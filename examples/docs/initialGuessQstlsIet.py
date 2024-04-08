import numpy as np
import pandas as pd
import qupled.quantum as qpq

# Define a QstlsIet object to solve a QSTLS-IET scheme
qstls = qpq.QstlsIet(30.0, 1.0, "QSTLS-HNC",
                     mixing = 0.2,
                     resolution = 0.1,
                     cutoff = 5,
                     matsubara = 16,
                     scheme2DIntegrals = "segregated",
                     threads = 16)

# Solve the scheme
qstls.compute()

# Create a custom initial guess from the output files of the previous run
qstls.setGuess("rs30.000_theta1.000_QSTLS-HNC.h5")

# Solve the scheme again with the new initial guess and coupling parameter
qstls.compute()
