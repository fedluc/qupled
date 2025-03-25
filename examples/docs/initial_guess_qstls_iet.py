from qupled.quantum import QstlsIet

# Define the object used to solve the scheme
qstls = QstlsIet()

# Define the input parameters
inputs = QstlsIet.Input(10.0, 1.0, "QSTLS-HNC")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.integral_strategy = "segregated"

# Solve the scheme
qstls.compute(inputs)

# Create a custom initial guess from the output files of the previous run
inputs.guess = QstlsIet.get_initial_guess("rs10.000_theta1.000_QSTLS-HNC.h5")

# Solve the scheme again with the new initial guess
qstls.compute(inputs)
