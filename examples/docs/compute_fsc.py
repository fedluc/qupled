import numpy as np
import qupled.stls as SCHEME   # or any other scheme
from qupled.output import DataBase
import qupled.fsc as fsc

N_VALUES   = [10, 100, 1000] # Number of particles
RS_TARGET  = 1.0
THETA      = 1.0
G_CUTOFF   = 50.0 # Should be kept above 50 for higher accuracy
DRS        = 0.1
SCHEME_KW  = {"mixing": 0.5, "cutoff": G_CUTOFF, "dimension": "D3"}  

inputs = fsc.Input(coupling=RS_TARGET, degeneracy=THETA, cutoff=G_CUTOFF, dimension=SCHEME_KW["dimension"], drs=DRS)
solver = fsc.Solver(scheme_module=SCHEME, scheme_kwargs=SCHEME_KW)
run_id = solver.compute(inputs, N_values=N_VALUES)

out = DataBase.read_results(run_id, names=["fsc_uint", "fsc_fxc"])
fsc_uint = np.atleast_1d(out.get("fsc_uint", np.array([])))
fsc_fxc  = np.atleast_1d(out.get("fsc_fxc",  np.array([])))
for N, u, f in zip(N_VALUES, fsc_uint, fsc_fxc):
    print(f"N={N:>6}  fsc_uint={u:.6e}  fsc_fxc={f:.6e}")