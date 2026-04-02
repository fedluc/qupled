import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from qupled.postprocess.output import DataBase, OutputType
from qupled.util.dimension import Dimension
from qupled.schemes import rpa, hf, stls

# Define an Rpa object to solve the RPA scheme
solver = stls.Solver()
solver.compute(stls.Input(1.0, 1.0, dimension=Dimension._3D))
solver.compute_itcf()

# Retrieve information from the output files
results = DataBase.read_results(OutputType.SCHEME, solver.run_id)

# Compare the results for the from the two schemes in a plot
colormap = cm["viridis"].reversed()
nl = results["itcf"].shape[1]
for i in range(nl):
    color = colormap(1.0 - 1.0 * i / nl)
    plt.plot(
        results["wvg"],
        results["itcf"][:, i],
        color=color,
    )
plt.show()
