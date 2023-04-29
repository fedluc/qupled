import numpy as np
import pandas as pd
import qupled.qupled as qp
import qupled.Static as Static

# Define a Qstls object to solve the QSTLS scheme
qstls = Static.Qstls(10.0, 1.0,
                     mixing = 0.5,
                     resolution = 0.1,
                     cutoff = 10,
                     matsubara = 16,
                     threads = 16)

# Solve the QSTLS scheme
qstls.compute()

# Create a custom initial guess from the output files of the previous run
guess = qp.QstlsGuess()
fileName = "rs10.000_theta1.000_QSTLS.h5"
guess.wvg = pd.read_hdf(fileName, "wvg")[0].to_numpy()
guess.ssf = pd.read_hdf(fileName, "ssf")[0].to_numpy()
qtls.qInputs.guess = guess

# Change the coupling parameter
qtls.input.coupling = 30.0

# Solve the scheme again with the new initial guess and coupling parameter
qtls.compute()
