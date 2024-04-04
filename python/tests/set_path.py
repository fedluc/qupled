import os
import sys

cwd = os.getcwd()
sys.path.insert(0, cwd)
sys.path.insert(1, os.path.join(cwd, "tests", "examples/docs"))
