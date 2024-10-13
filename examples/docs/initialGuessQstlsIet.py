from qupled.quantum import QstlsIet

# Define the input parameters
inputs = QstlsIet.Input(10.0, 1.0, "QSTLS-HNC")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.int2DScheme = "segregated"

# Solve the scheme
QstlsIet(inputs).compute()

# Create a custom initial guess from the output files of the previous run
guess = QstlsIet.getInitialGuess("rs10.000_theta1.000_QSTLS-HNC.h5")

# Solve the scheme again with the new initial guess
inputs.guess = guess
QstlsIet(inputs).compute()
