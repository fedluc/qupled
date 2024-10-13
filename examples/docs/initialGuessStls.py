from qupled.classic import Stls
from qupled.classic import IterativeScheme

# Define the input parameters
inputs = Stls.Input(10.0, 1.0)
inputs.mixing = 0.2

# Solve scheme
Stls(inputs).compute()

# Create a custom initial guess from the output files of the previous run
inputs.guess = Stls.getInitialGuess("rs10.000_theta1.000_STLS.h5")

# Solve the scheme again with the new initial guess and coupling parameter
Stls(inputs).compute()
