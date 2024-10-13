import numpy as np
import pandas as pd
from qupled.classic import Stls
from qupled.util import Hdf

# Define the input parameters
inputs = Stls.Input(10.0, 1.0)
inputs.mixing = 0.5

# Define an Stls object to solve the STLS scheme
stls = Stls(inputs)

# Solve scheme
stls.compute()

# Plot some results
stls.plot(["ssf", "slfc", "rdf"])

# Plot the ideal density response for a few matsubara frequencies
stls.plot(["idr"], matsubara=np.arange(1, 10, 2))

# Access the static structure factor from the Stls object
ssf = stls.ssf
print("Static structure factor from the Stls object: ")
print(ssf)

# Access the static structure factor from the output file
ssf = Hdf().read(stls.hdfFileName, ["ssf"])["ssf"]
print("Static structure factor from the output file: ")
print(ssf)

# Compute the internal energy
print("Internal energy: ")
print(stls.computeInternalEnergy())
