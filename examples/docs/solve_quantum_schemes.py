import qupled.qstls as qstls
import qupled.qstlsiet as qstlsiet

# Define a Qstls object to solve the QSTLS scheme
scheme = qstls.Qstls()

# Define the input parameters
inputs = qstls.Input(10.0, 1.0)
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16

# Solve the QSTLS scheme
scheme.compute(inputs)
print(scheme.results.uint)

# Define a QstlsIet object to solve the QSTLS-IET scheme
scheme = qstlsiet.QstlsIet()

# Define the input parameters for one of the QSTLS-IET schemes
inputs = qstlsiet.Input(10.0, 1.0, "QSTLS-LCT")
inputs.mixing = 0.5
inputs.matsubara = 16
inputs.threads = 16
inputs.integral_strategy = "segregated"

# solve the QSTLS-IET scheme
scheme.compute(inputs)
