import qupled.qstls as qstls

# Define the object used to solve the scheme
scheme = qstls.Qstls()

# Define the input parameters
inputs = qstls.Input(10.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16

# Solve the QSTLS scheme and store the internal energy (v1 calculation)
scheme.compute(inputs)
uInt1 = scheme.results.uint

# Pass in input the fixed component of the auxiliary density response
inputs.fixed = "adr_fixed_theta1.000_matsubara16_QSTLS.bin"

# Repeat the calculation and recompute the internal energy (v2 calculation)
scheme.compute(inputs)
uInt2 = scheme.results.uint

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)

# Change the coupling parameter
inputs.coupling = 20.0

# Compute with the updated coupling parameter
scheme.compute(inputs)

# Change the degeneracy parameter
inputs.degeneracy = 2.0

# Compute with the update degeneracy parameter (this throws a RunTimeError)
scheme.compute(inputs)
