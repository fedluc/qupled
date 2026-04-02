from qupled.schemes import qstls, qstlsiet

# Define a Qstls object to solve the QSTLS scheme
scheme = qstls.Solver()

# Define the input parameters
inputs = qstls.Input(10.0, 1.0, mixing=0.5, matsubara=16, threads=16)

# Solve the QSTLS scheme
scheme.compute(inputs)
