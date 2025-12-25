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

# Analyze outputs
for run_id in run_ids:
    # Fetch finite size correction results
    run_data = output.DataBase.read_run(
        output.OutputType.FINITE_SIZE_CORRECTION, run_id
    )
    run = run_data["run"]
    results = run_data["results"]
    inputs = run_data["inputs"]
    # Fetch scheme results
    scheme_result = output.DataBase.read_results(
        output.OutputType.SCHEME, run["scheme_run_id"]
    )
    # Print summary
    print("--------------------------------")
    print("Finite Size Correction Results:")
    print(f"Number of Particles: {inputs['number_of_particles']}")
    print(f"Corrections: uint = {results['uint']:.6e}, fxc = {results['fxc']:.6e}")
