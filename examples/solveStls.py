import numpy as np
import pandas as pd
import qupled.classic as qpc

# Define an Stls object to solve the STLS scheme
stls = qpc.Stls(10.0, # Coupling parameter
                1.0,  # Degeneracy parameter
                mixing = 0.5,
                cutoff = 10)

# Solve scheme
stls.compute()

# Plot some results
stls.plot(["ssf", "slfc", "rdf"])
