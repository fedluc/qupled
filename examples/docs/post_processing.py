import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from qupled.schemes import stls

# Solve the STLS scheme in the ground state for rs=10.0
solver = stls.Solver()
solver.compute(stls.Input(10.0, 0.0, mixing=0.7))

# Check the default values for the quantites to be computed in the post-processing step
print(f"----- Default values for the post-processing quantities -----")
print(f"Radial distribution function: {solver.results.rdf}")
print(f"Imaginary-time correlation function: {solver.results.itcf}")
print("-------------------------------------------------------------")

# Compute the radial distribution function and the imaginary-time correlation function
solver.compute_rdf()
solver.compute_itcf(tau=[0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])

# Plot the radial distribution function
plt.plot(solver.results.rdf_grid, solver.results.rdf, color="b")
plt.xlabel("Interparticle distance")
plt.ylabel("Radial distribution function")
plt.show()

# Plot the imaginary-time correlation function
colormap = cm["viridis"].reversed()
# nl = min(solver.results.itcf.shape[1], 10)
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
