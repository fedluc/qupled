from qupled.classic import StlsIet

# Define the input parameters
inputs = StlsIet.Input(10.0, 1.0, "STLS-HNC")
inputs.mixing = 0.5

# Solve scheme with HNC bridge function
StlsIet(inputs).compute()

# Change to a dielectric scheme with a different bridge function
inputs.theory = "STLS-LCT"

# Solve again with an LCT bridge function
StlsIet(inputs).compute()
