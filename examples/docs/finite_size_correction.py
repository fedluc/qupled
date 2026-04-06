from qupled.postprocess import output, finite_size_correction as fsc
from qupled.schemes import rpa

# Define the scheme
solver = rpa.Solver()
inputs = rpa.Input(coupling=1.0, degeneracy=1.0, cutoff=50.0)

# Compute the finite size correction for different number of particles
run_ids = []
finite_size_correction = fsc.FiniteSizeCorrection()
for number_of_particles in [10, 100, 1000]:
    fsc_inputs = fsc.Input(
        number_of_particles=number_of_particles,
        drs=0.1,
        scheme=inputs,
    )
    finite_size_correction.compute(solver, fsc_inputs)
    run_ids.append(finite_size_correction.run_id)

# Fetch and print the results
col1, col2, col3 = 21, 27, 26
header = f"│ {'Number of Particles':>{col1}} │ {'uint correction':>{col2}} │ {'fxc correction':>{col3}} │"
top = f"┌{'─' * (col1 + 2)}┬{'─' * (col2 + 2)}┬{'─' * (col3 + 2)}┐"
mid = f"├{'─' * (col1 + 2)}┼{'─' * (col2 + 2)}┼{'─' * (col3 + 2)}┤"
bot = f"└{'─' * (col1 + 2)}┴{'─' * (col2 + 2)}┴{'─' * (col3 + 2)}┘"
print("\nFinite Size Correction Results:")
print(top)
print(header)
print(mid)
for run_id in run_ids:
    run_data = output.DataBase.read_run(
        run_id, type=output.OutputType.FINITE_SIZE_CORRECTION
    )
    number_of_particles = run_data.inputs["number_of_particles"]
    uint = run_data.results["uint"]
    fxc = run_data.results["fxc"]
    print(f"│ {number_of_particles:>{col1}} │ {uint:>{col2}.6e} │ {fxc:>{col3}.6e} │")
print(bot)
