import numpy as np
import pandas as pd
import qupled.Static as Static

# Define an Stls object to solve the STLS scheme
stls = Static.Stls(10.0, # Coupling parameter
                   1.0,  # Degeneracy parameter
                   mixing = 0.5,
                   cutoff = 10)

# Solve scheme
stls.compute()

# Plot some results
stls.plot(["ssf", "slfc", "rdf"])

# Plot the ideal density response for a few matsubara frequencies
stls.plot(["idr"], matsubara = np.arange(1, 10, 2))

# Access the static structure factor from the Stls object
ssf = stls.scheme.ssf
print("Static structure factor from the Stls object: ")
print(ssf)

# Access the static structure factor from the output file
ssf =  pd.read_hdf("rs10.000_theta1.000_STLS.h5", "ssf")[0].to_numpy()
print("Static structure factor from the output file: ")
print(ssf)
