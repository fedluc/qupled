from pprint import pprint
import matplotlib.pyplot as plt
from qupled.classic import Rpa, ESA
import qupled.util as qpu

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme #########")
rpa = Rpa()

# Solve the RPA scheme
rpa.compute(Rpa.Input(10.0, 1.0))

# Define an ESA object to solve the ESA scheme
print("######### Solving the ESA scheme #########")
esa = ESA()

# Solve the ESA scheme
esa.compute(ESA.Input(10.0, 1.0))

# Inspect the outuput files to see what data was saved
outputFileRPA = rpa.hdf_file_name
outputFileESA = esa.hdf_file_name
print("########## Data stored for the RPA scheme #########")
pprint(qpu.HDF().inspect(outputFileRPA))
print("########## Data stored for the ESA scheme #########")
pprint(qpu.HDF().inspect(outputFileRPA))

# Retrieve some information that we want to plot from the output files
hdfDataRPA = qpu.HDF().read(
    outputFileRPA, ["coupling", "degeneracy", "theory", "ssf", "wvg"]
)
hdfDataESA = qpu.HDF().read(
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
