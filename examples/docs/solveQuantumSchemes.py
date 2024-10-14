import numpy as np
from qupled.quantum import Qstls, QstlsIet

# Define the input parameters
inputs = Qstls.Input(10.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16

# Define a Qstls object to solve the QSTLS scheme
qstls = Qstls(inputs)

# Solve the QSTLS scheme
qstls.compute()

# Plot the density responses and the static local field correction
qstls.plot(["idr", "adr", "slfc"], matsubara=np.arange(1, 10, 2))

# Define the input parameters for one of the QSTLS-IET schemes
inputs = QstlsIet.Input(10.0, 1.0, "QSTLS-LCT")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.int2DScheme = "segregated"

# Define a QstlsIet object to solve the QSTLS-IET scheme
qstls = QstlsIet(inputs)

# solve the QSTLS-IET scheme
qstls.compute()

