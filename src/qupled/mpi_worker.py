import json
from pathlib import Path


def run_mpi_worker(InputCls, ResultCls, NativeInputCls, NativeSchemeCls):
    with Path("input.json").open() as f:
        input_dict = json.load(f)
    inputs = InputCls.from_dict(input_dict)
    native_inputs = NativeInputCls()
    inputs.to_native(native_inputs)
    scheme = NativeSchemeCls(native_inputs)
    scheme.compute()
    results = ResultCls()
    results.from_native(scheme)
    with Path("results.json").open("w") as f:
        json.dump(results.to_dict(), f)
