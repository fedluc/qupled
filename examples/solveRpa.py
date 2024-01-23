import numpy as np
import pandas as pd
import qupled.classic as qpc

# Define an Rpa object to solve the RPA scheme
rpa = qpc.Rpa(10.0, # Coupling parameter
              1.0,  # Degeneracy parameter
              cutoff = 10)

# Solve scheme
rpa.compute()

# Plot some results
rpa.plot(["ssf", "slfc", "rdf"])

# Plot the ideal density response for a few matsubara frequencies
rpa.plot(["idr"], matsubara = np.arange(1, 10, 2))

# Access the static structure factor from the Stls object
ssf = rpa.scheme.ssf
print("Static structure factor from the Rpa object: ")
print(ssf)

# Access the static structure factor from the output file
ssf =  pd.read_hdf("rs10.000_theta1.000_RPA.h5", "ssf")[0].to_numpy()
print("Static structure factor from the output file: ")
print(ssf)
