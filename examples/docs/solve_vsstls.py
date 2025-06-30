import qupled.vsstls as vsstls

# Define the object used to solve the scheme
scheme = vsstls.Solver()

# Define the input parameters
inputs = vsstls.Input(2.0, 1.0, mixing=0.5, alpha=[-0.2, 0.2])

# Compute
scheme.compute(inputs)

# Load the free energy integrand computed for rs = 2.0
fxci = scheme.get_free_energy_integrand(scheme.run_id)

# Setup a new VSStls simulation for rs = 5.0
inputs.coupling = 5.0
inputs.alpha = [0.5, 0.7]
inputs.free_energy_integrand = fxci

# Compute
scheme.compute(inputs)
