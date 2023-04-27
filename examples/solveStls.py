import numpy as np
import qupled.Static as Static

# Define an Stls object to solve the STLS scheme
stls = Static.Stls(10.0, # Coupling parameter
                   1.0,  # Degeneracy parameter
                   mixing = 0.5,
                   cutoff = 10)

# Solve scheme
stls.compute()

# Plot some results
stls.plot(["ssf", "slfc", "rdf"])

# Plot the ideal density response for a few matsubara frequencies
stls.plot(["idr"], matsubara = np.arange(1, 10, 2))

# Access directly the computed static structure factor
ssf = stls.scheme.ssf
print(ssf)
