from qupled.classic import Stls

# Define the object used to solve the scheme
stls = Stls()

# Define the input parameters
inputs = Stls.Input(10.0, 1.0)
inputs.mixing = 0.2

# Solve scheme
stls.compute(inputs)

# Create a custom initial guess from the output files of the previous run
inputs.guess = Stls.get_initial_guess("rs10.000_theta1.000_STLS.h5")

# Solve the scheme again with the new initial guess and coupling parameter
stls.compute(inputs)
