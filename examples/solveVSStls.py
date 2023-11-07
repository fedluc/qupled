import pandas as pd
import qupled.qupled as qp
import qupled.classic as qpc

# Define a VSStls object to solve the VS-STLS scheme
stls = qpc.VSStls(10.0, 
                  0.0,
                  mixing = 0.8,
                  resolution = 0.1,
                  cutoff = 10,
                  couplingResolution = 0.1,
                  alpha = [0.5, 0.8],
                  iterations = 100,
                  threads = 9)

# Compute
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci"])

# Setup a new VSStls simulation for rs = 15.0 and use the free energy
# integrand computed for rs = 10.0
stls.inputs.coupling = 15.0
fxci = qp.FreeEnergyIntegrand()
fileName = "rs10.000_theta0.000_VSSTLS.h5"
fxci.grid = pd.read_hdf(fileName, "fxcGrid")[0].to_numpy()
fxci.integrand = pd.read_hdf(fileName, "fxci")[0].to_numpy()
stls.inputs.freeEnergyIntegrand = fxci

# Compute
stls.compute()

# Plot the results
stls.plot(["ssf", "slfc", "fxci"])

