import numpy as np
import pandas as pd
import qupled.classic as qpc

# Define an ESA object to solve the ESA scheme
ESA = qpc.ESA(15.0, # Coupling parameter
              1.0,  # Degeneracy parameter
              cutoff = 10)

# Solve scheme
ESA.compute()

# Plot some results
ESA.plot(["ssf", "slfc", "rdf"])


# Access the static structure factor from the output file
slfc =  pd.read_hdf("rs15.000_theta1.000_ESA.h5", "slfc")[0].to_numpy()
ssf = pd.read_hdf("rs15.000_theta1.000_ESA.h5", "ssf")[0].to_numpy()
rdf = pd.read_hdf("rs15.000_theta1.000_ESA.h5", "rdf")[0].to_numpy()
print("Static local field correction from the output file: ")
print(slfc)
print("Static structure factor from the output file: ")
print(ssf)
print("Radial distribution function from the output file: ")
print(rdf)
