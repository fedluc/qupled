import qupled.qvsstls as qvsstls

# Define the object used to solve the scheme
scheme = qvsstls.QVSStls()

# Define the input parameters
inputs = qvsstls.Input(1.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.alpha = [-0.2, 0.4]
inputs.iterations = 100
inputs.threads = 16

# Solve the QVSSTLS scheme and store the internal energy (v1 calculation)
scheme.compute(inputs)
uInt1 = scheme.results.uint

# Pass in input the fixed component of the auxiliary density response
inputs.fixed = "adr_fixed_theta1.000_matsubara16_QVSSTLS.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
scheme.compute(inputs)
uInt2 = scheme.results.uint

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
