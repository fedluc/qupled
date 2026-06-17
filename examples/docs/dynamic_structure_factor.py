import matplotlib.pyplot as plt
import numpy as np

from qupled.schemes import rpa

# Solve a finite-temperature, three-dimensional static scheme.
theta = 1.0
inputs = rpa.Input(
    coupling=2.0,
    degeneracy=theta,
    cutoff=3.0,
    resolution=0.2,
    frequency_cutoff=20.0,
)
solver = rpa.Solver()
solver.compute(inputs)

# Compute the dynamic structure factor on a positive-frequency grid.
frequency = np.arange(0.0, inputs.frequency_cutoff + 0.05, 0.05)
solver.compute_dsf(frequency)

# Compute the ITCF directly and by integrating the DSF.
beta = 1.0 / theta
tau = beta * np.linspace(0.0, 1.0, 11)
solver.compute_itcf(tau)
solver.compute_itcf_from_dsf(tau)

# Select the computed wave vector closest to the requested value.
target_x = 1.0
q_index = np.argmin(np.abs(solver.results.wvg - target_x))
x = solver.results.wvg[q_index]

# Plot S(x, Omega) for one wave vector.
plt.plot(solver.results.frequency, solver.results.dsf[q_index])
plt.xlabel(r"$\Omega$")
plt.ylabel(r"$S(x,\Omega)$")
plt.title(rf"$x={x:.1f}$")
plt.show()

# Compare the ITCF obtained directly with the ITCF reconstructed from the DSF.
plt.plot(
    solver.results.tau / beta,
    solver.results.itcf[q_index],
    label="Matsubara",
)
plt.plot(
    solver.results.tau / beta,
    solver.results.itcf_from_dsf[q_index],
    "--",
    label="DSF integral",
)
plt.xlabel(r"$\tau/\beta$")
plt.ylabel(r"$F(x,\tau)$")
plt.title(rf"$x={x:.1f}$")
plt.legend()
plt.show()
