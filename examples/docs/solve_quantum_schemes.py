import numpy as np
from qupled.quantum import Qstls, QstlsIet

# Define a Qstls object to solve the QSTLS scheme
qstls = Qstls()

# Define the input parameters
inputs = Qstls.Input(10.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16

# Solve the QSTLS scheme
qstls.compute(inputs)

# Plot the density responses and the static local field correction
qstls.plot(["idr", "adr", "slfc"], matsubara=np.arange(1, 10, 2))

# Define a QstlsIet object to solve the QSTLS-IET scheme
qstls = QstlsIet()

# Define the input parameters for one of the QSTLS-IET schemes
inputs = QstlsIet.Input(10.0, 1.0, "QSTLS-LCT")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.integral_strategy = "segregated"

# solve the QSTLS-IET scheme
qstls.compute(inputs)
