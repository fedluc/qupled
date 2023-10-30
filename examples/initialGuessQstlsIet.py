import numpy as np
import pandas as pd
import qupled.qupled as qp
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
guess = qp.QstlsGuess()
guess = qp.QstlsGuess()
fileName = "rs30.000_theta1.000_QSTLS-HNC.h5"
guess.wvg = pd.read_hdf(fileName, "wvg")[0].to_numpy()
guess.ssf = pd.read_hdf(fileName, "ssf")[0].to_numpy()
guess.adr = np.ascontiguousarray(pd.read_hdf(fileName, "adr").to_numpy())
guess.matsubara = pd.read_hdf(fileName, "inputs")["matsubara"][0].tolist()
qstls.qInputs.guess = guess

# Change the coupling parameter
qstls.inputs.coupling = 20.0

# Solve the scheme again with the new initial guess and coupling parameter
qstls.compute()
