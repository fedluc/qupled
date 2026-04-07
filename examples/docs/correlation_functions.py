import matplotlib.pyplot as plt
from qupled.util.dimension import Dimension
from qupled.postprocess.correlation_functions import compute_rdf
from qupled.postprocess.output import DataBase
from qupled.schemes import stls

# Solve an STLS scheme
scheme = stls.Solver()
scheme.compute(stls.Input(2.0, 1.0, dimension=Dimension._2D, mixing=0.7))
run_id = scheme.run_id

# After solving, the RDF is not yet in the database
run = DataBase.read_run(run_id, result_names=["rdf", "rdf_grid"])
print(
    "Radial distribution function available before computation: ", "rdf" in run.results
)

# Compute the RDF and store it back under the same run_id
compute_rdf(run_id)

# The RDF is now available in the database
run = DataBase.read_run(run_id)
print(
    "Radial distribution function available after the computation: ",
    "rdf" in run.results,
)
plt.plot(run.results["rdf_grid"], run.results["rdf"], color="b")
plt.xlabel("Interparticle distance")
plt.ylabel("Radial distribution function")
plt.show()
