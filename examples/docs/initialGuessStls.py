import numpy as np
import pandas as pd
import qupled.classic as qpc

# Define an Stls object to solve the STLS scheme
stls = qpc.Stls(10.0,
                1.0,
                mixing = 0.2,
                cutoff = 10)

# Solve scheme
stls.compute()

# Create a custom initial guess from the output files of the previous run
stls.setGuess("rs10.000_theta1.000_STLS.h5")

# Change the coupling parameter
stls.inputs.coupling = 10.0

# Solve the scheme again with the new initial guess and coupling parameter
stls.compute()
