import numpy as np
import qupled.quantum as qpq

# Define a QVSStls object to solve the QVSStls scheme
qvsstls = qpq.QVSStls(
    1.0,
    1.0,
    mixing=0.5,
    matsubara=16,
    couplingResolution=0.1,
    degeneracyResolution=0.1,
    alpha=[-0.2, 0.4],
    errorIntegrals=1e-5,
    iterations=100,
    threads=16,
)

# Solve the QVSSTLS scheme and store the internal energy (v1 calculation)
qvsstls.compute()
uInt1 = qvsstls.computeInternalEnergy()

# Pass in input the fixed component of the auxiliary density response
qvsstls.inputs.fixed = "adr_fixed_theta1.000_matsubara16.zip"

# Repeat the calculation and recompute the internal energy (v2 calculation)
qvsstls.compute()
uInt2 = qvsstls.computeInternalEnergy()

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)
