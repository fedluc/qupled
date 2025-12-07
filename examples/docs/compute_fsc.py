from qupled import rpa
from qupled.output import DataBase
import qupled.finite_size_correction as fsc

finite_size_correction = fsc.FiniteSizeCorrection()
solver = rpa.Solver()
inputs = fsc.Input(
    scheme = rpa.Input(coupling=1.0, degeneracy=1.0, cutoff=50.0),
    drs=0.1,
    number_of_particles=100,
)
run_id = finite_size_correction.compute(solver, inputs)
print(run_id)
# solver = fsc.Solver(scheme_module=SCHEME, scheme_kwargs=SCHEME_KW)
# run_id = solver.compute(inputs, N_values=N_VALUES)

# out = DataBase.read_results(run_id, names=["fsc_uint", "fsc_fxc"])
# fsc_uint = np.atleast_1d(out.get("fsc_uint", np.array([])))
# fsc_fxc = np.atleast_1d(out.get("fsc_fxc", np.array([])))
# for N, u, f in zip(N_VALUES, fsc_uint, fsc_fxc):
#     print(f"N={N:>6}  fsc_uint={u:.6e}  fsc_fxc={f:.6e}")
