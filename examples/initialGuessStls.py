import numpy as np
import pandas as pd
import qupled.qupled as qp
import qupled.Static as Static

# Define an Stls object to solve the STLS scheme
stls = Static.Stls(10.0,
                   1.0,
                   mixing = 0.2,
                   cutoff = 10)

# Solve scheme
stls.compute()

# Create a custom initial guess from the output files of the previous run
guess = qp.SlfcGuess()
fileName = "rs10.000_theta1.000_STLS.h5"
guess.wvg = pd.read_hdf(fileName, "wvg")[0].to_numpy()
guess.slfc = pd.read_hdf(fileName, "slfc")[0].to_numpy()
stls.inputs.guess = guess

# Change the coupling parameter
stls.inputs.coupling = 30.0

# Solve the scheme again with the new initial guess and coupling parameter
stls.compute()
