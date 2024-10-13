import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from qupled.classic import Rpa, ESA
import qupled.util as qpu

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme #########")
rpa = Rpa(Rpa.Input(10.0, 1.0))

# Solve the RPA scheme
rpa.compute()

# Define an ESA object to solve the ESA scheme
print("######### Solving the ESA scheme #########")
esa = ESA(ESA.Input(10.0, 1.0))

# Solve the ESA scheme
esa.compute()

# Inspect the outuput files to see what data was saved
outputFileRPA = rpa.hdfFileName
outputFileESA = esa.hdfFileName
print("########## Data stored for the RPA scheme #########")
print(yaml.dump(qpu.Hdf().inspect(outputFileRPA)))
print("########## Data stored for the ESA scheme #########")
print(yaml.dump(qpu.Hdf().inspect(outputFileESA)))

# Retrieve some information that we want to plot from the output files
hdfDataRPA = qpu.Hdf().read(
    outputFileRPA, ["coupling", "degeneracy", "theory", "ssf", "wvg"]
)
hdfDataESA = qpu.Hdf().read(
    outputFileESA, ["coupling", "degeneracy", "theory", "ssf", "wvg"]
)

# Compare the results for the from the two schemes in a plot
plt.plot(hdfDataRPA["wvg"], hdfDataRPA["ssf"], color="b", label=hdfDataRPA["theory"])
plt.plot(hdfDataESA["wvg"], hdfDataESA["ssf"], color="r", label=hdfDataESA["theory"])
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title(
    "State point : (coupling = "
    + str(hdfDataRPA["coupling"])
    + ", degeneracy = "
    + str(hdfDataRPA["degeneracy"])
    + ")"
)
plt.show()
