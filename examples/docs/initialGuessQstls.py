import numpy as np
import pandas as pd
import qupled.quantum as qpq

# Define a Qstls object to solve the QSTLS scheme
qstls = qpq.Qstls(
    10.0, 1.0, mixing=0.4, resolution=0.1, cutoff=10, matsubara=16, threads=16
)

# Solve the QSTLS scheme
qstls.compute()

# Create a custom initial guess from the output files of the previous run
qstls.setGuess("rs10.000_theta1.000_QSTLS.h5")

# Change the coupling parameter
qstls.inputs.coupling = 10.0

# Solve the scheme again with the new initial guess and coupling parameter
qstls.compute()
