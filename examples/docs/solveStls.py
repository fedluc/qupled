import numpy as np
import pandas as pd
import qupled.classic as qpc
import qupled.util as qpu

# Define an Stls object to solve the STLS scheme
stls = qpc.Stls(10.0, 1.0, mixing=0.5, cutoff=10)

# Solve scheme
stls.compute()

# Plot some results
stls.plot(["ssf", "slfc", "rdf"])

# Plot the ideal density response for a few matsubara frequencies
stls.plot(["idr"], matsubara=np.arange(1, 10, 2))

# Access the static structure factor from the Stls object
ssf = stls.scheme.ssf
print("Static structure factor from the Stls object: ")
print(ssf)

# Access the static structure factor from the output file
ssf = qpu.Hdf().read(stls.hdfFileName, ["ssf"])["ssf"]
print("Static structure factor from the output file: ")
print(ssf)

print(stls.computeInternalEnergy())
