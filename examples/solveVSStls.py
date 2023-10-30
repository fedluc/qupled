import pandas as pd
import qupled.qupled as qp
import qupled.classic as qpc

# Define a VSStls object to solve the VS-STLS scheme
stls = qpc.VSStls(0.2, 
                  0.0,
                  mixing = 0.8,
                  resolution = 0.1,
                  cutoff = 5,
                  couplingResolution = 0.01,
                  alpha = 0.0,
                  mixingAlpha = 1.0)

# Compute
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci"])

# Setup a new VSStls simulation  and use the free energy
# integrand computed for rs = 0.2
fxci = qp.FreeEnergyIntegrand()
fileName = "rs0.200_theta0.000_VSSTLS.h5"
fxci.grid = pd.read_hdf(fileName, "fxcGrid")[0].to_numpy()
fxci.integrand = pd.read_hdf(fileName, "fxci")[0].to_numpy()
stls.inputs.freeEnergyIntegrand = fxci

# Compute up to rs = 0.4
stls.inputs.coupling = 0.4
stls.inputs.mixingAlpha = 0.7
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci"])
