import matplotlib.pyplot as plt
from qupled.postprocess.output import DataBase
from qupled.schemes import esa, rpa

# Solve two schemes: RPA and ESA
scheme = rpa.Solver()
scheme.compute(rpa.Input(10.0, 1.0))
run_id_rpa = scheme.run_id

scheme = esa.Solver()
scheme.compute(esa.Input(10.0, 1.0))
run_id_esa = scheme.run_id

# Retrieve information from the database
rpa_data = DataBase.read_run(run_id_rpa)
esa_data = DataBase.read_run(run_id_esa)
rpa_results = rpa_data["results"]
rpa_inputs = rpa_data["inputs"]
esa_results = esa_data["results"]
esa_inputs = esa_data["inputs"]

# Plot the static structure factor for the two schemes
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
