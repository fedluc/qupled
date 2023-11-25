import pandas as pd
import numpy as np
import qupled.qupled as qp
import qupled.classic as qpc

# Define a VSStls object to solve the VS-STLS scheme
stls = qpc.VSStls(5.0, 
                  1.0,
                  mixing = 0.5,
                  resolution = 0.1,
                  cutoff = 10,
                  couplingResolution = 0.1,
                  degeneracyResolution = 0.01,
                  alpha = [-0.2, 0.2],
                  errorIntegrals = 1e-5,
                  iterations = 100,
                  threads = 9)


# Compute
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci", "sdr"])

# Setup a new VSStls simulation for rs = 10.0 and use the free energy
# integrand computed for rs = 5.0
stls.inputs.coupling = 10.0
stls.inputs.alpha = [0.5, 0.7]
stls.setFreeEnergyIntegrand("rs5.000_theta1.000_VSSTLS.h5")

# Compute
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci"])
