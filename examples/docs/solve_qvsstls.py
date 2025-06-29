import qupled.qvsstls as qvsstls

# Define the object used to solve the scheme
scheme = qvsstls.Solver()

# Define the input parameters
inputs = qvsstls.Input(1.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.alpha = [-0.2, 0.4]
inputs.iterations = 100
inputs.threads = 16

# Solve scheme for rs = 1.0
scheme.compute(inputs)

# Load the free energy integrand computed for rs = 1.0
fxci = scheme.get_free_energy_integrand(scheme.run_id)

# Setup a new  simulation for rs=2.0
inputs.coupling = 2.0
inputs.alpha = [0.1, 0.5]
inputs.free_energy_integrand = fxci

# Solve scheme for rs = 2.0
scheme.compute(inputs)
