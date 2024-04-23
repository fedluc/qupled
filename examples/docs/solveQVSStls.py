import numpy as np
import qupled.quantum as qpq

# Solve the QVSStls scheme for rs=1.0 and theta=1.0
qvsstls = qpq.QVSStls(1.0, 
                      1.0,
                      mixing = 0.5,
                      matsubara = 16,
                      alpha = [-0.2, 0.4],
                      errorIntegrals = 1e-5,
                      iterations = 100,
                      threads = 16)

# Solve the QVSSTLS scheme
qvsstls.compute()

# Setup a new  simulation for rs=2.0 and theta=1.0 and use the free energy
# integrand computed for rs = 1.0
qvsstls.inputs.coupling = 2.0
qvsstls.inputs.alpha = [0.1, 0.5]
qvsstls.setFreeEnergyIntegrand("rs1.000_theta1.000_QVSSTLS.h5")
qvsstls.compute()

# Plot the results
qvsstls.plot(["ssf", "fxci"])

