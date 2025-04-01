import qupled.stls as stls
from qupled.output import DataBase

# Define the object used to solve the scheme
scheme = stls.Stls()

# Define the input parameters
inputs = stls.Input(10.0, 1.0)
inputs.mixing = 0.5

# Solve scheme
scheme.compute(inputs)

# Access the internal energy from the output file
results = DataBase().read_results(scheme.run_id, names=["uint"])
print("Internal energy from the output file: ")
print(results["uint"])

# Compute the internal energy
print("Internal energy from the result class: ")
print(scheme.results.uint)
