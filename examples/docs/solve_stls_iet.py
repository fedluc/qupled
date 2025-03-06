from qupled.classic import StlsIet

# Define the object used to solve the scheme
stls = StlsIet()

# Define the input parameters
inputs = StlsIet.Input(10.0, 1.0, "STLS-HNC")
inputs.mixing = 0.5

# Solve scheme with HNC bridge function
stls.compute(inputs)

# Change to a dielectric scheme with a different bridge function
inputs.theory = "STLS-LCT"

# Solve again with an LCT bridge function
stls.compute(inputs)
