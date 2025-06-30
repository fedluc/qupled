import qupled.stls as stls

# Define the object used to solve the scheme
scheme = stls.Solver()

# Define the input parameters
inputs = stls.Input(10.0, 1.0, mixing=0.2)

# Solve scheme
scheme.compute(inputs)

# Create a custom initial guess from the output files of the previous run
inputs.guess = scheme.get_initial_guess(scheme.run_id)

# Solve the scheme again with the new initial guess and coupling parameter
scheme.compute(inputs)
