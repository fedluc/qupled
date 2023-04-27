import qupled.Static as Static

# Define an StlsIet object to solve an STLS-IET scheme
stls = Static.StlsIet(10.0, 
                      1.0,
                      "STLS-HNC",
                      mixing = 0.5,
                      cutoff = 10)

# Solve scheme with HNC bridge function
stls.compute()

# Change to a dielectric scheme with a different bridge function
stls.inputs.theory = "STLS-LCT"

# Solve again
stls.compute()
