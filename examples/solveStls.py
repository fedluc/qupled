import qupled.qupled as qp
import qupled.Static as Static

# Define an Stls object to solve STLS scheme
stls = Static.Stls(10.0, # Coupling parameter
                   1.0,  # Degeneracy parameter
                   mixing = 0.5,
                   cutoff = 5)

# Solve scheme
stls.compute()

# Plot solutions
stls.plot(["ssf", "sdr", "slfc", "ssfHF", "bf", "idr", "rdf"])

