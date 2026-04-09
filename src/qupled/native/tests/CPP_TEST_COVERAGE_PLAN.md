# Native C++ Test Coverage Plan

## Scope I inspected

I reviewed all C++ sources and headers under `src/qupled/native`:

- `33` production `.cpp` files
- `36` `.hpp` files
- `14` current test `.cpp` files

This plan focuses on full coverage of the native codebase in that tree.

## File inventory and test status

### 1) Utility layer (`include/util`, `src/util`)

Production files:

- `src/qupled/native/src/util/database.cpp`
- `src/qupled/native/src/util/dimensions_util.cpp`
- `src/qupled/native/src/util/logger.cpp`
- `src/qupled/native/src/util/mpi_util.cpp`
- `src/qupled/native/src/util/num_util.cpp`
- `src/qupled/native/src/util/numerics.cpp`
- `src/qupled/native/src/util/vector2D.cpp`
- `src/qupled/native/src/util/vector3D.cpp`
- `src/qupled/native/src/util/vector_util.cpp`

Coverage direction:

- Keep and expand low-level deterministic tests (vector math, tolerance helpers, dimensional dispatch).
- Add edge/error-path tests for interpolators/integrators/solvers (`numerics`).
- Add DB adapter tests around `database.cpp` behavior and expected failure modes.

### 2) Thermodynamic layer (`include/thermo`, `src/thermo`)

Production files:

- `src/qupled/native/src/thermo/chemical_potential.cpp`
- `src/qupled/native/src/thermo/free_energy.cpp`
- `src/qupled/native/src/thermo/internal_energy.cpp`
- `src/qupled/native/src/thermo/itcf.cpp`
- `src/qupled/native/src/thermo/rdf.cpp`
- `src/qupled/native/src/thermo/thermo_util.cpp`

Coverage direction:

- Keep unit tests for free/internal energy, RDF, and chemical potential.
- Add branch-focused tests for `itcf.cpp` and `thermo_util.cpp`:
  - empty-input guards
  - HF early-return branch
  - finite vs ground-state branches
  - size-mismatch error paths

### 3) Scheme core (`include/schemes`, `src/schemes`)

Production files:

- `src/qupled/native/src/schemes/esa.cpp`
- `src/qupled/native/src/schemes/hf.cpp`
- `src/qupled/native/src/schemes/iet.cpp`
- `src/qupled/native/src/schemes/input.cpp`
- `src/qupled/native/src/schemes/qstls.cpp`
- `src/qupled/native/src/schemes/qstlsiet.cpp`
- `src/qupled/native/src/schemes/qvsstls.cpp`
- `src/qupled/native/src/schemes/rpa.cpp`
- `src/qupled/native/src/schemes/stls.cpp`
- `src/qupled/native/src/schemes/stlsiet.cpp`
- `src/qupled/native/src/schemes/vsstls.cpp`

Coverage direction:

- `input.cpp`: mostly covered; add remaining corner cases for default/sentinel flows.
- `hf.cpp`/`rpa.cpp`: add compact-grid regression tests for finite and ground-state dispatch paths.
- `stls.cpp`/`stlsiet.cpp`: convergence-step and initial-guess path tests.
- `esa.cpp`: coefficient caching and invariant checks (finite outputs, monotonic/physical bounds where applicable).
- `iet.cpp`: bridge function mapping/theory routing + domain errors.
- `qstls*` and `qvs*`: staged tests with temp sqlite/blob fixtures and small grids to verify:
  - read/write roundtrips
  - fixed-component loading branches
  - matrix-shape consistency and no-crash on minimal valid inputs.

### 4) VS infrastructure (`include/vs`, `src/vs`)

Production files:

- `src/qupled/native/src/vs/vsbase.cpp`
- `src/qupled/native/src/vs/vsmanager.cpp`

Coverage direction:

- Add targeted tests for finite-difference stencils and derivative metadata mapping.
- Verify update/application ordering invariants in manager loops.
- Use lightweight fake workers to unit-test `VSManager` logic without heavy physics kernels.

### 5) Python interface bindings (`include/python_interface`, `src/python_interface`)

Production files:

- `src/qupled/native/src/python_interface/inputs.cpp`
- `src/qupled/native/src/python_interface/native.cpp`
- `src/qupled/native/src/python_interface/schemes.cpp`
- `src/qupled/native/src/python_interface/util.cpp`
- `src/qupled/native/src/python_interface/utilities.cpp`

Coverage direction:

- Add native-side conversion tests for `python_interface/util.cpp` (array shape/order checks, conversion correctness).
- Add binding smoke tests to verify module/class/function registration is alive.
- If needed, split these into a dedicated binding test target to keep `native-cpp` unit runtime stable.

## Execution plan (phased)

## Phase 1.5 - Test tree reorganization (completed)

- Decision: mirror test folders to production module layout before expanding Phase 2+ coverage.
- Implemented structure:
  - `tests/common` (global env / suite bootstrap)
  - `tests/util`
  - `tests/thermo`
  - `tests/schemes`
  - `tests/vs`
  - `tests/fixtures`
- Build wiring updated so tests can include shared helpers via `../tests` include path.
- Rationale: reduce churn for upcoming additions and make per-module coverage gaps obvious.

## Phase 0 - Test infrastructure hardening

- Keep global test environment for MPI init/finalize.
- Add reusable fixtures/helpers:
  - minimal valid `Input`/`StlsInput`/`QstlsInput` builders
  - temporary sqlite/blob fixture utilities
  - shared numeric tolerance helpers

## Phase 1 - Finish deterministic foundations

- Complete utility + thermo branch coverage first (fast and stable).
- Ensure all pure math/integration tests are deterministic and tolerance-bounded.

## Phase 2 - Scheme behavior matrix

- Add per-scheme tests in increasing complexity:
  - `HF` -> `RPA` -> `STLS` -> `STLS-IET` -> `ESA`
- For each scheme test:
  - success path on minimal valid input
  - at least one key error path
  - one invariant assertion on output shape/value sanity

## Phase 3 - Quantum and persistence-heavy schemes

- Add isolated tests for `Qstls`, `QstlsIet`, `QVSStls` with temp db/files.
- Verify database interaction branches and roundtrip correctness.
- Keep computational grids tiny to avoid slow/fragile tests.

## Phase 4 - VS orchestration and bindings

- Add `VSManager`/`VSBase` unit tests with fake workers.
- Add pybind smoke/conversion tests.

## Phase 5 - Coverage gate and cleanup

- Run full suite via `./devtool build --native-only --native-tests` and `./devtool test native-cpp`.
- Remove overlaps/duplication, standardize naming and fixture usage.
- Optionally add coverage reporting target if desired.

## Proposed test folder structure

- `src/qupled/native/tests/util/*_test.cpp`
- `src/qupled/native/tests/thermo/*_test.cpp`
- `src/qupled/native/tests/schemes/*_test.cpp`
- `src/qupled/native/tests/vs/*_test.cpp`
- `src/qupled/native/tests/python_interface/*_test.cpp`
- `src/qupled/native/tests/fixtures/*` (helpers, builders, temp-db utilities)

This keeps tests discoverable by subsystem and makes ownership/maintenance clearer as the suite grows.
