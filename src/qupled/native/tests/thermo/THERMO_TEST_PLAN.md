# Thermo Module Test Plan

## Goal

Build extensive, maintainable coverage for all classes in `src/qupled/native/include/thermo`, plus high-level free functions in `thermoUtil`.

This plan follows the same structure and discipline used for schemes tests and is intended to be reused in future sessions.

## Scope

Headers in scope:

- `thermo/free_energy.hpp`
- `thermo/internal_energy.hpp`
- `thermo/rdf.hpp`
- `thermo/chemical_potential.hpp`
- `thermo/itcf.hpp`
- `thermo/thermo_util.hpp`

Source companions in `src/qupled/native/src/thermo/*.cpp` are covered via public-interface tests and util-level behavior tests.

## Testing principles

- Test every public method of every public class.
- Test both top-level classes and `thermoUtil` free functions as first-class interfaces.
- Keep tests deterministic: tiny grids, explicit tolerances, no hidden randomness.
- Assert both value behavior and shape/invariant contracts.
- Validate error paths as rigorously as success paths.

## Test granularity strategy

- Default to one test per behavior (pytest-style), even in C++.
- Prefer one test per public method branch (success, invalid input, boundary case).
- Keep each test focused on a single failure reason for faster diagnosis.
- Group assertions only when they represent one inseparable behavior and share expensive setup.
- Treat repeated construction of different objects in the same test as a refactor signal: split into independent tests unless the setup/assertion is inseparable.
- Prefer separate tests for each dimension/branch split (2D vs 3D, `tau=0` vs `tau>0`, `x=0` vs `x>0`, coupling zero vs finite coupling).

## Constructor coverage strategy

- Add dedicated constructor tests whenever constructors initialize critical state or feed branch-sensitive behavior.
- For classes with guarded construction logic (if any), include explicit `RejectsInvalid...` and `AcceptsValid...` tests.
- Keep constructor tests independent from heavier compute-path tests.

## What to do with private methods

- Do **not** test private methods directly.
- Cover private logic through public observable behavior.
- If complex private logic becomes hard to reach, extract testable helpers instead of using visibility hacks.

## Suggested test file layout

Under `src/qupled/native/tests/thermo/`:

- `free_energy_test.cpp`
- `internal_energy_test.cpp`
- `rdf_test.cpp`
- `chemical_potential_test.cpp`
- `itcf_test.cpp`
- `thermo_util_test.cpp`

## Detailed coverage matrix

## 1) `free_energy.hpp`

Class:

- `FreeEnergy`

Coverage:

- `get()` for non-normalized integration (dedicated test).
- `get()` for normalized integration (dedicated test).
- `rs == 0` + normalized branch (`-inf`).

## 2) `internal_energy.hpp`

Class:

- `InternalEnergy`

Coverage:

- `get()` in 2D and 3D (separate tests).
- Zero signal for `S(k)=1` baseline.
- Deterministic non-zero case sanity against known closed-form setup.

## 3) `rdf.hpp`

Class:

- `Rdf`

Coverage:

- `get()` in 2D and 3D via separate tests.
- `r=0` and `r>0` branches via separate tests.
- `S(k)=1` baseline (`g(r)=1`).

## 4) `chemical_potential.hpp`

Class:

- `ChemicalPotential`

Coverage:

- 2D closed-form branch.
- 3D root-solver branch with residual check.
- pre-compute sentinel behavior (`get()` before compute).

## 5) `itcf.hpp`

Classes:

- `thermoUtil::ItcfNonInteracting`
- `thermoUtil::ItcfNonInteractingGround`
- `thermoUtil::Itcf`
- `thermoUtil::ItcfGround`

Coverage:

- Non-interacting finite: `tau=0` and `tau>0` branches.
- Non-interacting ground: `tau=0`, `x=0`, and standard branch behavior, each in dedicated tests.
- Interacting finite (`Itcf`): `x=0`, zero-coupling fallback, `tau=0`, and `tau>0` paths in dedicated tests.
- Interacting ground (`ItcfGround`): `x=0`, zero-coupling fallback, `tau=0`, and `tau>0` paths in dedicated tests.

## 6) `thermo_util.hpp`

Functions:

- `computeInternalEnergy`
- `computeFreeEnergy` (both overloads)
- `computeRdf`
- `computeItcf`

Coverage:

- Return-shape/value sanity for each function.
- Overload behavior (`normalize=true/false`) for free energy.
- out-of-range free-energy error path.
- ITCF empty-input short-circuit.
- ITCF HF branch bypassing LFC validation (dedicated test).
- ITCF interacting-branch LFC size validation error (dedicated test).

## Execution order

1. Scalar observables (`FreeEnergy`, `InternalEnergy`, `Rdf`, `ChemicalPotential`).
2. ITCF classes (`ItcfNonInteracting`, `ItcfNonInteractingGround`, `Itcf`, `ItcfGround`).
3. `thermoUtil` free-function integration paths.
4. Final sweep for missing public methods and error paths.

## Quality gates for this module

- Every public method in `include/thermo/*.hpp` has at least one direct test.
- Every documented error branch has at least one failure-path test.
- `./devtool build --native-only --native-tests` and `./devtool test native-cpp` pass.

## Execution environment note

- Run validation commands inside the devcontainer, not on the host machine.
- Current container: `4deb4a17d562`.
- Example:
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool build --native-only --native-tests'`
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool test native-cpp'`

## Commit author policy

- Author for all commits touching this thermo test plan/execution must be:
  - `Federico Lucco Castello <71817492+fedluc@users.noreply.github.com>`
