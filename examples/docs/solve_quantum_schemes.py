import qupled.qstls as qstls
import qupled.qstlsiet as qstlsiet

# Define a Qstls object to solve the QSTLS scheme
scheme = qstls.Solver()

# Define the input parameters
inputs = qstls.Input(10.0, 1.0, mixing=0.5, matsubara=16, threads=16)

# Solve the QSTLS scheme
scheme.compute(inputs)
print(scheme.results.uint)

# Define a QstlsIet object to solve the QSTLS-IET scheme
scheme = qstlsiet.Solver()

# Define the input parameters for one of the QSTLS-IET schemes
inputs = qstlsiet.Input(
    10.0,
    1.0,
    theory="QSTLS-LCT",
    mixing=0.5,
    matsubara=16,
    threads=16,
    integral_strategy="segregated",
)

# solve the QSTLS-IET scheme
scheme.compute(inputs)
