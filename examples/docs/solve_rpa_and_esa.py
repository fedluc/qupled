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
rpa_run = DataBase.read_run(run_id_rpa)
esa_run = DataBase.read_run(run_id_esa)


# Plot the static structure factor for the two schemes
plt.plot(
    rpa_run.results["wvg"],
    rpa_run.results["ssf"],
    color="b",
    label=rpa_run.inputs["theory"],
)
plt.plot(
    esa_run.results["wvg"],
    esa_run.results["ssf"],
    color="r",
    label=esa_run.inputs["theory"],
)
plt.legend(loc="lower right")
plt.xlabel("Wave vector")
plt.ylabel("Static structure factor")
plt.title(
    "State point : (coupling = "
    + str(rpa_run.inputs["coupling"])
    + ", degeneracy = "
    + str(rpa_run.inputs["degeneracy"])
    + ")"
)
plt.show()
