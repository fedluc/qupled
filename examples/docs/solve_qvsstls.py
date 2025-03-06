from qupled.quantum import QVSStls

# Define the object used to solve the scheme
qvsstls = QVSStls()

# Define the input parameters
inputs = QVSStls.Input(1.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.alpha = [-0.2, 0.4]
inputs.iterations = 100
inputs.threads = 16

# Solve scheme for rs = 1.0
qvsstls.compute(inputs)

# Load the free energy integrand computed for rs = 1.0
fxci = QVSStls.get_free_energy_integrand("rs1.000_theta1.000_QVSSTLS.h5")

# Setup a new  simulation for rs=2.0
inputs.coupling = 2.0
inputs.alpha = [0.1, 0.5]
inputs.free_energy_integrand = fxci

# Solve scheme for rs = 2.0
qvsstls.compute(inputs)

# Plot the results
qvsstls.plot(["ssf", "fxc_int"])
