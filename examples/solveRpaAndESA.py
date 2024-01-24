import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qupled.classic as qpc

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme")
rpa = qpc.Rpa(10.0, # Coupling parameter
              1.0,  # Degeneracy parameter
              cutoff = 10)

# Solve the RPA scheme
rpa.compute()

# Define an ESA object to solve the ESA scheme
print("######### Solving the ESA scheme")
esa = qpc.ESA(10.0,
              1.0,
              cutoff = 10)

# Solve the ESA scheme
esa.compute()

# Retrieve the input information for the two schemes
outputFileRPA = "rs10.000_theta1.000_RPA.h5"
outputFileESA = "rs10.000_theta1.000_ESA.h5"
inputRPA = pd.read_hdf(outputFileRPA, "inputs")
inputESA = pd.read_hdf(outputFileESA, "inputs")
rs = inputRPA["coupling"].to_string(index=False)
theta = inputRPA["degeneracy"].to_string(index=False)
theoryRPA = inputRPA["theory"].to_string(index=False)
theoryESA = inputESA["theory"].to_string(index=False)

# Access the wave-vector grid
wvg = pd.read_hdf(outputFileRPA, "wvg")[0].to_numpy()

# Access the static structure factor
ssfRpa =  pd.read_hdf(outputFileRPA, "ssf")[0].to_numpy()
ssfESA = pd.read_hdf(outputFileESA, "ssf")[0].to_numpy()

# Compare the results from the two schemes
plt.plot(wvg, ssfRpa, color="b")
plt.plot(wvg, ssfESA, color="r")
plt.legend([theoryRPA, theoryESA], loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title("State point : (coupling = " + rs + ", degeneracy = " + theta + ")")
plt.show()
