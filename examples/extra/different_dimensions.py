import matplotlib.pyplot as plt
from qupled.postprocess.output import DataBase, OutputType
from qupled.util.dimension import Dimension
from qupled.schemes import rpa

# Define an Rpa object to solve the RPA scheme
print("######### Solving the RPA scheme in 3D #########")
rpa3D_scheme = rpa.Solver()
rpa3D_scheme.compute(rpa.Input(10.0, 1.0, dimension=Dimension._3D))

# Define an ESA object to solve the ESA scheme
print("######### Solving the RPA scheme in 2D #########")
rpa2D_scheme = rpa.Solver()
rpa2D_scheme.compute(rpa.Input(10.0, 1.0, dimension=Dimension._2D))

# Retrieve information from the output files
run_3D = DataBase.read_run(rpa3D_scheme.run_id)
run_2D = DataBase.read_run(rpa2D_scheme.run_id)

# Compare the results for the from the two schemes in a plot
plt.plot(
    run_3D.results["wvg"],
    run_3D.results["ssf"],
    color="b",
    label=Dimension.from_dict(run_3D.inputs["dimension"]).value,
)
plt.plot(
    run_2D.results["wvg"],
    run_2D.results["ssf"],
    color="r",
    label=Dimension.from_dict(run_2D.inputs["dimension"]).value,
)
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title(
    "State point : (coupling = "
    + str(run_3D.inputs["coupling"])
    + ", degeneracy = "
    + str(run_3D.inputs["degeneracy"])
    + ")"
)
plt.show()
