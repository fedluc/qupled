from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
from qupled.schemes import hf
from qupled.util.dimension import Dimension

# Solve two schemes: RPA and ESA
solver = hf.Solver()
solver.compute(hf.Input(1.0, 0.0, dimension=Dimension._2D))

# Plot the imaginary-time correlation function
colormap = cm["viridis"].reversed()
nl = min(solver.results.idr.shape[1], 10)
# nl = solver.results.idr.shape[1]
for i in range(nl):
    color = colormap(1.0 - 1.0 * i / nl)
    plt.plot(
        solver.results.wvg,
        solver.results.idr[:, i],
        color=color,
        label=f"$l={i}$",
    )
plt.legend(loc="center right")
plt.xlabel("Wave-vector")
plt.ylabel("Density response function")
plt.savefig("test_2D_ground.png")
