from qupled.classic import Stls
from qupled.util import DataBase

# Define the object used to solve the scheme
stls = Stls()

# Define the input parameters
inputs = Stls.Input(10.0, 1.0)
inputs.mixing = 0.5

# Solve scheme
stls.compute(inputs)

# Access the internal energy from the output file
results = DataBase().read_results(stls.run_id, names=["uint"])
print("Internal energy from the output file: ")
print(results["uint"])

# Compute the internal energy
print("Internal energy from the result class: ")
print(stls.results.uint)
