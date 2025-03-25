from qupled.quantum import QstlsIet

# Define the object used to solve the scheme
qstls = QstlsIet()

# Define the input parameters
inputs = QstlsIet.Input(10.0, 1.0, "QSTLS-HNC")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.integral_strategy = "segregated"

# Solve the QSTLS-HNC scheme and store the internal energy (v1 calculation)
qstls.compute(inputs)
uInt1 = qstls.compute_internal_energy()

# Pass in input the fixed component of the auxiliary density response
inputs.fixed = "adr_fixed_theta1.000_matsubara16_QSTLS-HNC.bin"
inputs.fixed_iet = "adr_fixed_theta1.000_matsubara16_QSTLS-HNC.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
qstls.compute(inputs)
uInt2 = qstls.compute_internal_energy()

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
