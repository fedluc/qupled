import qupled.stlsiet as stlsiet

# Define the object used to solve the scheme
scheme = stlsiet.Solver()

# Define the input parameters
inputs = stlsiet.Input(10.0, 1.0, theory="STLS-HNC")
inputs.mixing = 0.5

# Solve scheme with HNC bridge function
scheme.compute(inputs)

# Change to a dielectric scheme with a different bridge function
inputs.theory = "STLS-LCT"

# Solve again with an LCT bridge function
scheme.compute(inputs)
