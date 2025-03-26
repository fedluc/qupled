import matplotlib.pyplot as plt
from qupled.classic import Rpa, ESA
from qupled.util import DataBase

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme #########")
rpa = Rpa()
rpa.compute(Rpa.Input(10.0, 1.0))

# Define an ESA object to solve the ESA scheme
print("######### Solving the ESA scheme #########")
esa = ESA()
esa.compute(ESA.Input(10.0, 1.0))

# Retrieve information from the output files
rpa_data = DataBase.read_run(rpa.run_id)
esa_data = DataBase.read_run(esa.run_id)
rpa_results = rpa_data["results"]
rpa_inputs = rpa_data["inputs"]
esa_results = esa_data["results"]
esa_inputs = esa_data["inputs"]

# Compare the results for the from the two schemes in a plot
plt.plot(rpa_results["wvg"], rpa_results["ssf"], color="b", label=rpa_inputs["theory"])
plt.plot(esa_results["wvg"], esa_results["ssf"], color="r", label=esa_inputs["theory"])
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title(
    "State point : (coupling = "
    + str(rpa_inputs["coupling"])
    + ", degeneracy = "
    + str(rpa_inputs["degeneracy"])
    + ")"
)
plt.show()
