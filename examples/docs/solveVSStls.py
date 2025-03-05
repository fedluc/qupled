from qupled.classic import VSStls

# Define the object used to solve the scheme
vsstls = VSStls()

# Define the input parameters
inputs = VSStls.Input(5.0, 1.0)
inputs.mixing = 0.5
inputs.degeneracy_resolution = 0.01
inputs.alpha = [-0.2, 0.2]

# Compute
vsstls.compute(inputs)

# Plot the results
vsstls.plot(["ssf", "slfc", "fxci", "sdr"])

# Load the free energy integrand computed for rs = 5.0
fxci = VSStls.get_free_energy_integrand("rs5.000_theta1.000_VSSTLS.h5")

# Setup a new VSStls simulation for rs = 10.0
inputs.coupling = 10.0
inputs.alpha = [0.5, 0.7]
inputs.free_energy_integrand = fxci

# Compute
vsstls.compute(inputs)

# Plot the results
vsstls.plot(["ssf", "slfc", "fxci"])
