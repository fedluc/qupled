# MPI Implementation Analysis

This note reviews how MPI is orchestrated in qupled. It separates the Python
process orchestration from the native C++ implementation because they have
different failure modes and different public API implications.

## Python MPI Orchestration

### Current Design

The public user-facing interface is the solver API:

```python
inputs = stls.Input(coupling=10.0, degeneracy=1.0, processes=2)
solver = stls.Solver()
solver.compute(inputs)
```

The MPI subprocess protocol is private plumbing. `hf.Solver._compute_native_mpi`
writes `input.json`, launches `mpiexec -n N python -m <scheme-module>`, then
reads `status.json` and `results.json`. Each scheme module has an
`if __name__ == "__main__"` block that calls `Solver.run_mpi_worker(...)` with
the appropriate input and result classes.

The authored docs and examples recommend setting `processes` on the input
dataclass. They do not recommend direct use of commands such as
`python -m qupled.schemes.hf`.

### Strengths

- The public interface is simple: users only set `processes > 1` when the
  native extension is built with MPI support and `mpiexec` is available.
- Running MPI in subprocesses gives each rank a fresh Python process and keeps
  MPI lifecycle concerns out of the already-running parent Python process.
- The root-rank-only result write is straightforward and avoids merging Python
  objects across ranks.

### Weaknesses

- ~~MPI is selected solely by `native.uses_mpi`, not by `inputs.processes`.
  `hf.Solver._compute_native()` currently takes the MPI path whenever the native
  extension was built with MPI, even when `processes == 1`. The docstring says
  it should require more than one process.~~
- The subprocess command uses the literal `"python"` instead of
  `sys.executable`, which can invoke the wrong interpreter in virtualenvs,
  editable installs, CI, or HPC modules.
- The file protocol uses fixed names in the current working directory:
  `input.json`, `results.json`, and `status.json`. Concurrent runs, stale files,
  interrupted runs, or nested solver calls can collide.
- Cleanup happens only after all reads succeed. A failed subprocess or failed
  result read leaves the fixed temporary files behind.
- Result and status are written as separate non-atomic JSON files. A crash can
  leave partial output or mismatched status/results.
- Input reconstruction uses `serializable_dataclass.from_dict`, which constructs
  objects with `cls.__new__(cls)`. That bypasses dataclass initialization,
  default factories, and `__post_init__` validation in MPI workers.
- VS-type schemes mutate shared input objects during subcalls and restore state
  only on success. A subprocess failure can leave `inputs.coupling` or related
  precomputation fields mutated.

### Recommendations

1. Keep the public user interface unchanged: `Input(..., processes=N)` followed
   by `solver.compute(inputs)`.
2. ~~Select MPI only when both conditions are true:
   `native.uses_mpi and inputs.processes > 1`.~~
3. Replace the per-scheme `python -m qupled.schemes.<scheme>` worker entry
   points with one centralized private worker module, for example
   `python -m qupled.util.mpi_worker --solver qupled.schemes.hf:Solver ...`.
   Scheme-level `if __name__ == "__main__"` blocks can remain temporarily as
   compatibility shims.
4. Use `sys.executable` for subprocesses.
5. Use a per-run temporary directory, pass paths explicitly via CLI args or
   environment variables, and clean it up with `try/finally`.
6. Write one atomic output file, or write result/status files through temporary
   files followed by `Path.replace()`.
7. Rework dataclass deserialization so it calls `cls(**converted_values)` after
   recursive conversion. That preserves default handling and validation.
8. Use `try/finally` around VS precomputation state mutation.

## Native C++ MPI and Error Handling

### Current Design

The native layer centralizes MPI details in `MPIUtil`:

- `MPIUtil::init()` and `MPIUtil::finalize()` call `MPI_Init` and
  `MPI_Finalize` when `USE_MPI` is enabled.
- `MPIUtil::rank()`, `numberOfRanks()`, `isRoot()`, and `barrier()` wrap
  `MPI_COMM_WORLD`.
- `MPIUtil::parallelFor()` splits a loop across MPI ranks and optionally uses
  OpenMP inside each rank.
- `MPIUtil::gatherLoopData()` gathers rank-local loop output back into the full
  buffer on all ranks.
- `MPIUtil::throwError()` is used throughout the native code for validation,
  numerical failures, database/file failures, and unsupported feature paths.

The pybind wrapper calls `MPIUtil::init()` at the start of every native
`compute()` call and `MPIUtil::finalize()` at the end. Most solver `compute()`
methods catch `std::runtime_error`, print the message, and return `1`; success
returns `0`.

### Strengths

- The native code has a single MPI abstraction, so most scheme code does not
  include direct MPI calls.
- Root-only logging through `Logger` keeps normal progress output readable in
  MPI runs.
- The loop partitioning and gather pattern is simple and fits the existing
  wave-vector-grid computations.
- GSL errors are mostly routed through wrapper helpers instead of relying on
  the default GSL abort behavior.

### Main Error-Handling Weaknesses

#### 1. `MPIUtil::throwError()` mixes two different responsibilities

`MPIUtil::throwError()` throws a `std::runtime_error` only in single-process
mode. In multi-rank MPI mode it prints the message and calls `MPI_Abort`.

That makes the same validation error behave differently depending on process
count. Examples include invalid input values, unsupported dimensions, numerical
integration failures, missing files, and database errors.

This is the central weakness. Many errors are ordinary user/input/runtime
failures that should become a structured failed status or a Python exception.
Instead, under multi-rank MPI they bypass the solver `compute()` catch blocks,
prevent `status.json` from being written, and leave the parent Python process
with only a subprocess failure and stderr.

Recommended direction:

- Introduce a normal exception type, e.g. `QupledError`, for validation,
  numerical, and file/database errors.
- Keep MPI job abortion separate, e.g. `MPIUtil::abortJob(message)` or
  `MPIUtil::failCollectively(message)`.
- Use normal exceptions before collective computation starts.
- Use collective failure handling only when continuing would deadlock other
  ranks or corrupt a collective operation.

#### 2. Solver `compute()` catches too narrowly and converts errors too early

`HF::compute()` and `VSBase::compute()` catch only `std::runtime_error`, print
to stderr, and return `1`.

This has several consequences:

- `std::logic_error`, `std::out_of_range`, `std::bad_alloc`, and other
  `std::exception` subclasses can bypass the status path.
- Error details are lost at the Python API boundary except for stderr.
- Internal callers can ignore the returned status and continue with partial
  state.

There are already examples of ignored statuses:

- `Rpa::computeSsfHF()` calls `hf.compute()` and then consumes `hf.getSsf()`
  without checking the return value.
- `VSBase::alphaDifference()` calls `runGrid()` and ignores the returned
  status.

Recommended direction:

- Let internal solver code throw exceptions.
- Convert exceptions to status codes only at the outermost pybind/Python worker
  boundary.
- If status codes remain inside C++, require callers to check them and propagate
  failures immediately.
- Catch `const std::exception &` at the boundary, not just
  `std::runtime_error`.
- Preserve an error message alongside the numeric status.

#### 3. MPI lifecycle is not ownership-safe

`MPIUtil::init()` unconditionally calls `MPI_Init`, and `MPIUtil::finalize()`
unconditionally calls `MPI_Finalize`. The pybind wrapper runs this pair around
every `compute()`.

Risks:

- Calling `compute()` twice in the same process after MPI has been finalized is
  not safe.
- If another library, `mpi4py`, or a host application initialized MPI first,
  qupled may finalize MPI that it does not own.
- If an uncaught exception escapes between `init()` and `finalize()`,
  finalization is skipped.
- Hybrid OpenMP/MPI runs do not request or check MPI thread support.

Recommended direction:

- Track MPI ownership: call `MPI_Initialized` and `MPI_Finalized` before
  initializing/finalizing.
- Only finalize MPI if qupled initialized it.
- Use an RAII guard around compute-time MPI lifecycle.
- Consider `MPI_Init_thread` and verify the returned thread support level if MPI
  calls can occur from OpenMP regions.

#### 4. OpenMP exceptions are unsafe

`MPIUtil::parallelFor()` executes `loopFunc(i)` inside an OpenMP `parallel for`.
If `loopFunc` throws in an OpenMP worker thread, the behavior is not a normal
catchable exception path. In practice this can terminate the process. In
multi-rank mode, `throwError()` can also call `MPI_Abort` from an OpenMP worker
thread, without confirmed MPI thread support.

Recommended direction:

- Catch exceptions inside the OpenMP loop body.
- Store the first error message/status in thread-safe local state.
- After the parallel region, rethrow on the rank's main thread or enter a
  collective failure path.
- Avoid direct MPI calls from OpenMP worker threads unless MPI thread support is
  explicitly requested and checked.

#### 5. `gatherLoopData()` validates too late and ignores MPI return codes

`gatherLoopData()` calls `MPI_Allgatherv` before checking whether
`dataToGather` is null. It also relies on MPI's default error handling rather
than checking return codes.

Recommended direction:

- Validate `dataToGather`, `countsPerLoop`, and `loopData` before any MPI call.
- Use the same communicator wrapper consistently instead of mixing
  `MPICommunicator` and `MPI_COMM_WORLD`.
- Consider setting `MPI_ERRORS_RETURN` and wrapping MPI return codes in a helper
  that includes the MPI error string.
- Give null-buffer failures a useful error message.

#### 6. Constructor and input-setting errors are outside the status path

Several user-facing failures can happen before solver `compute()` starts:

- Native input setters validate values and call `throwError()`.
- Scheme constructors reject unsupported configurations such as unsupported 2D
  quantum schemes or ground-state QVSSTLS.
- `HF` builds the wave-vector grid in the constructor.

In the Python serial path, these become Python exceptions. In the MPI worker
path, construction happens in each worker process before the solver's
`compute()` catch block can run. If a constructor fails, no structured status or
result file is written.

Recommended direction:

- Validate inputs in Python before launching MPI where possible.
- Wrap the entire worker body, including native input conversion and native
  scheme construction, in one error boundary.
- Write a structured failure output from rank 0 when construction fails.

#### 7. Assertions are used for runtime assumptions

Several paths use `assert()` for assumptions that can be affected by user input,
state from files, or numerical edge cases. Examples include vector shape
assumptions and finite-difference bounds.

Assertions are not a good user-facing error mechanism: they abort the process
when enabled and disappear when `NDEBUG` is defined.

Recommended direction:

- Keep `assert()` for impossible internal invariants.
- Use explicit checks with useful exceptions for input/file/numerical state.

#### 8. Missing precomputed data can silently produce default data

`Qstls::readAdrFixed()` returns without error when the database query has no
matching row. The caller then continues with the caller-provided `Vector3D`
contents, which are often zero-initialized.

Recommended direction:

- Treat missing fixed-ADR rows as an explicit error that includes `runId` and
  ADR name.
- Validate binary file size before accepting the data.
- Make the database/file error path consistent in single-rank and multi-rank
  runs.

### Recommended Refactor Shape

The cleanest design is to make exceptions the native internal error mechanism
and reserve integer statuses for the Python boundary.

Suggested structure:

1. Add `QupledError`, possibly with categories such as `InputError`,
   `NumericsError`, `DatabaseError`, and `MpiError`.
2. Replace general uses of `MPIUtil::throwError()` with normal exceptions.
3. Add explicit MPI helpers:
   - `MPIUtil::check(int mpiStatus, std::string context)`
   - `MPIUtil::abortJob(std::string message)`
   - `MPIUtil::allRanksOk(bool localOk)` or similar collective status helpers
4. Wrap pybind `compute()` in an RAII MPI guard and a broad
   `catch (const std::exception&)` boundary.
5. Return or expose both status code and error message to Python.
6. Update the Python MPI worker to always write a structured failure output
   from rank 0 when possible.

This preserves the current public Python workflow while making failures much
easier to diagnose and less dependent on process count.
