from qupled.quantum import QVSStls

# Define the input parameters
inputs = QVSStls.Input(1.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.alpha = [-0.2, 0.4]
inputs.iterations = 100
inputs.threads = 16

# Define a QVSStls object to solve the QVSStls scheme
qvsstls = QVSStls(inputs)

# Solve the QVSSTLS scheme and store the internal energy (v1 calculation)
qvsstls.compute()
uInt1 = qvsstls.computeInternalEnergy()

# Pass in input the fixed component of the auxiliary density response
inputs.fixed = "adr_fixed_theta1.000_matsubara16.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
qvsstls = QVSStls(inputs)
qvsstls.compute()
uInt2 = qvsstls.computeInternalEnergy()

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
