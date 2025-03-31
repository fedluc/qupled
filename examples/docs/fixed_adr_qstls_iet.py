import qupled.qstlsiet as qstlsiet

# Define the object used to solve the scheme
scheme = qstlsiet.QstlsIet()

# Define the input parameters
inputs = qstlsiet.Input(10.0, 1.0, "QSTLS-HNC")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.integral_strategy = "segregated"

# Solve the QSTLS-HNC scheme and store the internal energy (v1 calculation)
scheme.compute(inputs)
uInt1 = scheme.results.uint

# Pass in input the fixed component of the auxiliary density response
inputs.fixed = "adr_fixed_theta1.000_matsubara16_QSTLS-HNC.bin"
inputs.fixed_iet = "adr_fixed_theta1.000_matsubara16_QSTLS-HNC.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
scheme.compute(inputs)
uInt2 = scheme.results.uint

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
