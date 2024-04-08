import numpy as np
import qupled.quantum as qpq

# Define a QVSStls object to solve the QVSSTLS scheme
qvsstls = qpq.QVSStls(0.1,
                      1.0,
                      mixing = 0.5,
                      resolution = 0.1,
                      cutoff = 10,
                      matsubara = 32,
                      threads = 8)

# Solves the QVSSTLS scheme for rs=0.1, theta=1.0 
# If a fixed file isn't specified this will produce
# 3 files containing the adrFixed components
# theta - dtheta/ theta/ theta + dtheta
qvsstls.compute()

# Start a QVSStls simulation for rs=1.0 and theta=1.0
qvsstls = qpq.QVSStls(rs = 1.0, 
                      theta = 1.0,
                      mixing = 0.5,
                      resolution = 0.1,
                      cutoff = 10,
                      matsubara = 32,
                      couplingResolution = 0.1,
                      degeneracyResolution = 0.1,
                      alpha = [-0.2, 0.4],
                      errorIntegrals = 1e-5,
                      iterations = 100,
                      errorAlpha = 1e-3,
                      threads = 8,
                      fixed="adr_fixed_theta1.000_matsubara32.bin")

# Solve the QVSSTLS scheme
qvsstls.compute()

# Setup a new QVSStls simulation for rs=2.0 and theta=1.0
# Additionally use the previous free energy for faster computation
qvsstls.inputs.coupling = 2.0
qvsstls.inputs.alpha = [0.1, 0.5]
qvsstls.setFreeEnergyIntegrand("rs1.000_theta1.000_QVSSTLS.h5")
qvsstls.compute()

# Plot the results
qvsstls.plot(["ssf"])

