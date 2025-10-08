import matplotlib.pyplot as plt
import qupled.rpa as rpa
import qupled.esa as esa
from qupled.output import DataBase
from qupled.dimension import Dimension

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme in 3D #########")
rpa3D_scheme = rpa.Solver()
rpa3D_scheme.compute(rpa.Input(10.0, 1.0, dimension=Dimension._3D))

# Define an ESA object to solve the ESA scheme
print("######### Solving the RPA scheme in 2D #########")
rpa2D_scheme = rpa.Solver()
rpa2D_scheme.compute(rpa.Input(10.0, 1.0, dimension=Dimension._2D))

# Retrieve information from the output files
rpa3D_data = DataBase.read_run(rpa3D_scheme.run_id)
rpa2D_data = DataBase.read_run(rpa2D_scheme.run_id)
rpa3D_results = rpa3D_data["results"]
rpa3D_inputs = rpa3D_data["inputs"]
rpa2D_results = rpa2D_data["results"]
rpa2D_inputs = rpa2D_data["inputs"]

# Compare the results for the from the two schemes in a plot
plt.plot(
    rpa3D_results["wvg"],
    rpa3D_results["ssf"],
    color="b",
    label=Dimension.from_dict(rpa3D_inputs["dimension"]).value,
)
plt.plot(
    rpa2D_results["wvg"],
    rpa2D_results["ssf"],
    color="r",
    label=Dimension.from_dict(rpa2D_inputs["dimension"]).value,
)
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title(
    "State point : (coupling = "
    + str(rpa3D_inputs["coupling"])
    + ", degeneracy = "
    + str(rpa3D_inputs["degeneracy"])
    + ")"
)
plt.show()
