import qupled.quantum as qpq

# Define a Qstls object to solve the QSTLS-HNC scheme
qstls = qpq.QstlsIet(
    30.0,
    1.0,
    "QSTLS-HNC",
    mixing=0.2,
    resolution=0.1,
    cutoff=5,
    matsubara=16,
    scheme2DIntegrals="segregated",
    threads=16,
)

# Solve the QSTLS-HNC scheme and store the internal energy (v1 calculation)
qstls.compute()
uInt1 = qstls.computeInternalEnergy()

# Pass in input the fixed component of the auxiliary density response
qstls.inputs.fixed = "adr_fixed_theta1.000_matsubara16.bin"
qstls.inputs.fixediet = "adr_fixed_theta1.000_matsubara16_QSTLS-HNC.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
qstls.compute()
uInt2 = qstls.computeInternalEnergy()

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
