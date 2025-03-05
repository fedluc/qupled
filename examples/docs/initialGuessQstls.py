from qupled.quantum import Qstls

# Define the object used to solve the scheme
qstls = Qstls()

# Define the input parameters
inputs = Qstls.Input(10.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16

# Solve the QSTLS scheme
qstls.compute(inputs)

# Create a custom initial guess from the output files of the previous run
inputs.guess = Qstls.get_initial_guess("rs10.000_theta1.000_QSTLS.h5")

# Solve the scheme again with the new initial guess
qstls.compute(inputs)
