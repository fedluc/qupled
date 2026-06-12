from qupled.postprocess.output import DataBase, OutputType
from qupled.schemes import stls
from qupled.util.dimension import Dimension

# Solve the same scheme with different numbers of processes to demonstrate multiple MPI executions in the same Python session.
run_ids = []
for nproc in [1, 2, 4]:
    scheme = stls.Solver()
    inputs = stls.Input(10.0, 1.0, mixing=0.2, processes=nproc, dimension=Dimension._3D)
    scheme.compute(inputs)
    run_ids.append((nproc, scheme.run_id))

# Compare results
for nproc, run_id in run_ids:
    results = DataBase.read_results(run_id, names=["uint"])
    print(f"Processes: {nproc}, Internal energy: {results['uint']}")
