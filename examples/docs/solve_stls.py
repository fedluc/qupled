import qupled.stls as stls
from qupled.output import DataBase, OutputType

# Define the object used to solve the scheme
scheme = stls.Solver()

# Define the input parameters
inputs = stls.Input(10.0, 1.0, mixing=0.5)

# Solve scheme
scheme.compute(inputs)

# Access the internal energy from the database
result_type = OutputType.SCHEME
results = DataBase.read_results(result_type, scheme.run_id, names=["uint"])
print("Internal energy from the output file: ")
print(results["uint"])

# Access the internal energy from the result class
print("Internal energy from the result class: ")
print(scheme.results.uint)
