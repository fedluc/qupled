from qupled.classic import VSStls

# Define the object used to solve the scheme
vsstls = VSStls()

# Define the input parameters
inputs = VSStls.Input(2.0, 1.0)
inputs.mixing = 0.5
inputs.alpha = [-0.2, 0.2]

# Compute
vsstls.compute(inputs)

# Load the free energy integrand computed for rs = 2.0
fxci = VSStls.get_free_energy_integrand(vsstls.run_id)

# Setup a new VSStls simulation for rs = 5.0
inputs.coupling = 5.0
inputs.alpha = [0.5, 0.7]
inputs.free_energy_integrand = fxci

# Compute
vsstls.compute(inputs)
