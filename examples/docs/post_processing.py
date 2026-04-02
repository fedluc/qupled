import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from qupled.postprocess.output import DataBase, OutputType
from qupled.util.dimension import Dimension
from qupled.schemes import rpa, hf, stls

# Solve the STLS scheme for a 2D system with rs=5.0
solver = stls.Solver()
solver.compute(stls.Input(5.0, 1.0, dimension=Dimension._2D, mixing=0.5))

# Check the default values for the quantites to be computed in the post-processing step
print(f"----- Default values for the post-processing quantities -----")
print(f"Radial distribution function: {solver.results.rdf}")
print(f"Imaginary-time correlation function: {solver.results.itcf}")
print("-------------------------------------------------------------")

# Compute the radial distribution function and the imaginary-time correlation function
solver.compute_rdf()
solver.compute_itcf()

# Plot the radial distribution function
plt.plot(solver.results.rdf_grid, solver.results.rdf, color="b")
plt.xlabel("Interparticle distance")
plt.ylabel("Radial distribution function")
plt.show()

# Plot the imaginary-time correlation function
colormap = cm["viridis"].reversed()
nl = solver.results.itcf.shape[1]
for i in range(nl):
    color = colormap(1.0 - 1.0 * i / nl)
    plt.plot(
        solver.results.wvg,
        solver.results.itcf[:, i],
        color=color,
        label=f"$\\tau^*={solver.results.tau[i]:.2f}$",
    )
plt.legend(loc="center right")
plt.xlabel("Wave-vector")
plt.ylabel("Imaginary-time correlation function")
plt.show()
