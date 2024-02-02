import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import qupled.classic as qpc

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme #########")
rpa = qpc.Rpa(10.0, # Coupling parameter
              1.0,  # Degeneracy parameter
              cutoff = 10)

# Solve the RPA scheme
rpa.compute()

# Define an ESA object to solve the ESA scheme
print("######### Solving the ESA scheme #########")
esa = qpc.ESA(10.0,
              1.0,
              cutoff = 10)

# Solve the ESA scheme
esa.compute()

# Inspect the outuput files to see what data was saved
outputFileRPA = rpa.hdfFileName
outputFileESA = esa.hdfFileName
print("########## Data stored for the RPA scheme #########")
print(yaml.dump(qpc.Hdf().inspect(outputFileRPA)))
print("########## Data stored for the ESA scheme #########")
print(yaml.dump(qpc.Hdf().inspect(outputFileESA)))

# Retrieve some information that we want to plot from the output files
hdfDataRPA = qpc.Hdf().read(outputFileRPA, ["coupling", "degeneracy", "theory", "ssf", "wvg"])
hdfDataESA = qpc.Hdf().read(outputFileESA, ["coupling", "degeneracy", "theory", "ssf", "wvg"])

# Compare the results for the from the two schemes in a plot
plt.plot(hdfDataRPA["wvg"], hdfDataRPA["ssf"], color="b", label = hdfDataRPA["theory"])
plt.plot(hdfDataESA["wvg"], hdfDataESA["ssf"], color="r", label = hdfDataESA["theory"])
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title("State point : (coupling = " + str(hdfDataRPA["coupling"])
          + ", degeneracy = " + str(hdfDataRPA["degeneracy"]) + ")")
plt.show()
