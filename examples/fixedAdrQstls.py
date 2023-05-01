import qupled.Static as Static

# Define a Qstls object to solve the QSTLS scheme
qstls = Static.Qstls(10.0, 1.0,
                     mixing = 0.5,
                     resolution = 0.1,
                     cutoff = 10,
                     matsubara = 16,
                     threads = 16)

# Solve the QSTLS scheme and store the internal energy (v1 calculation)
qstls.compute()
uInt1 = qstls.computeInternalEnergy()

# Pass in input the fixed component of the auxiliary density response
qstls.qInputs.fixed = "adr_fixed_rs10.000_theta1.000_QSTLS.bin"

# Repeat the calculation and recompute the internal energy (v2 calculation)
qstls.compute()
uInt2 = qstls.computeInternalEnergy()

# Compare the internal energies obtained with the two methods
print("Internal energy (v1) = %.8f" % uInt1)
print("Internal energy (v2) = %.8f" % uInt2)

# Change the coupling parameter
qstls.inputs.coupling = 20.0

# Compute with the updated coupling parameter
qstls.compute()

# Change the degeneracy parameter
qstls.inputs.degeneracy = 2.0

# Compute with the update degeneracy parameter (this throws an error)
qstls.compute() 
