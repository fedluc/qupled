# Schemes Module Test Plan

## Goal

Build extensive, maintainable coverage for all classes in `src/qupled/native/include/schemes`, including all related utility namespaces (`HFUtil`, `RpaUtil`, `StlsUtil`, `StlsIetUtil`, `QstlsUtil`, `QstlsIetUtil`, `ESAUtil`, `IetUtil`).

This plan is intentionally module-focused and intended to be reused in future sessions.

## Scope

Headers in scope:

- `schemes/input.hpp`
- `schemes/hf.hpp`
- `schemes/rpa.hpp`
- `schemes/stls.hpp`
- `schemes/stlsiet.hpp`
- `schemes/iet.hpp`
- `schemes/esa.hpp`
- `schemes/qstls.hpp`
- `schemes/qstlsiet.hpp`
- `schemes/vsstls.hpp`
- `schemes/qvsstls.hpp`

Source companions in `src/qupled/native/src/schemes/*.cpp` are covered via public API tests and utility-object tests.

## Testing principles

- Test every public method of every public class.
- Test utility namespace classes (`*Util`) as first-class units, not as incidental side effects.
- Keep tests deterministic: tiny grids, explicit tolerances, no hidden randomness.
- Assert both correctness and shape/invariant contracts.
- Validate error paths as rigorously as success paths.

## Test granularity strategy

- Default to one test per behavior (pytest-style), even in C++.
- Prefer one test per public method branch (success, invalid input, boundary case).
- Keep each test focused on a single failure reason for faster diagnosis.
- Group assertions only when they describe one inseparable behavior and share expensive setup.
- For expensive solver setup, keep granular assertions but reuse local fixtures/helpers to avoid duplication.

## Constructor coverage strategy

- Add dedicated constructor tests whenever constructors enforce constraints, allocate key state, or initialize required invariants.
- For each such class, include at least one `RejectsInvalid...` and one `AcceptsValid...` constructor test.
- Add constructor-initialization tests for observable post-construction state (e.g., array/grid sizes) when meaningful.
- Keep constructor tests independent from `compute()` tests so construction failures are diagnosed directly.

## What to do with private methods

Recommended policy:

- Do **not** test private methods directly.
- Cover private logic through public/protected observable behavior.
- If private logic is complex and difficult to drive via public API, extract it into a small helper in a testable namespace (`detail` or util helper) and test that helper directly.
- Avoid hacks like `#define private public`.
- Use friend-based testing only as a rare last resort for critical legacy code where extraction is impractical.

Reasoning:

- Private methods are implementation details and should be refactorable without breaking tests.
- Public/API-level tests better encode behavior guarantees.

## Suggested test file layout

Under `src/qupled/native/tests/schemes/`:

- `input_api_test.cpp`
- `hf_api_and_hfutil_test.cpp`
- `rpa_api_and_rpautil_test.cpp`
- `stls_api_and_stlsutil_test.cpp`
- `stlsiet_api_and_stlsietutil_test.cpp`
- `iet_api_and_ietutil_test.cpp`
- `esa_api_and_esautil_test.cpp`
- `qstls_api_and_qstlsutil_test.cpp`
- `qstlsiet_api_and_qstlsietutil_test.cpp`
- `vsstls_api_test.cpp`
- `qvsstls_api_test.cpp`

## Detailed coverage matrix

## 1) `input.hpp`

Classes:

- `Input`
- `IterationInput`
- `IetInput`
- `QuantumInput`
- `VSInput`
- concrete combos (`StlsInput`, `StlsIetInput`, `QstlsInput`, `QstlsIetInput`, `VSStlsInput`, `QVSStlsInput`)

Coverage:

- Every setter success + getter roundtrip.
- Every setter invalid-value exception path.
- Theory/mapping whitelist validation.
- Guess and free-energy-integrand shape/consistency checks.
- Sentinel/default behavior where relevant.

## 2) `hf.hpp` + `HFUtil`

Classes:

- `HF`
- `HFUtil::Idr`
- `HFUtil::IdrGround`
- `HFUtil::Ssf`
- `HFUtil::SsfGround`

Coverage:

- `HF::compute()` in finite and ground-state modes.
- Getter contracts (`wvg`, `idf/idr`, `ssf`, `lfc`, `sdr`, `uint`, `chemical_potential`).
- 2D and 3D branch behavior.
- Utility class `get()` methods for edge points (`x=0`, `l=0`, boundary frequencies).
- Known analytic reference points for `SsfGround` and `IdrGround` where available.

## 3) `rpa.hpp` + `RpaUtil`

Classes:

- `Rpa`
- `RpaUtil::SsfBase`
- `RpaUtil::Ssf`
- `RpaUtil::SsfGround`

Coverage:

- `Rpa::init`, `compute`, finite/ground branches.
- Zero-coupling fallback to HF behavior.
- Utility SSF classes on representative inputs (finite outputs, expected limits).
- Matsubara summation stability and static/dynamic LFC branch behavior.

## 4) `stls.hpp` + `StlsUtil`

Classes:

- `Stls`
- `StlsUtil::SlfcBase`
- `StlsUtil::Slfc`
- `StlsUtil::dynamic_pointer_cast` helpers

Coverage:

- STLS iterative loop: initial guess, convergence update, error computation.
- `initialGuessFromInput` success/failure branches.
- Mixing parameter effects.
- `Slfc` 2D/3D behavior, special points (`x=0`, `y=0`, `x=y`).
- Dynamic cast helper success and failure exception paths.

## 5) `iet.hpp` + `IetUtil`

Classes:

- `Iet`
- `IetUtil::BridgeFunction`

Coverage:

- `Iet::init`, bridge-function computation, initial-guess interpolation branch.
- Bridge function theory dispatch (`HNC`, `IOI`, `LCT`) and invalid-theory error.
- Mapping variants (`sqrt`, `linear`, standard) + invalid/ground-state forbidden branch.
- IOI/LCT validity-range exceptions and representative finite-value checks.

## 6) `stlsiet.hpp` + `StlsIetUtil`

Classes:

- `StlsIet`
- `StlsIetUtil::Slfc`

Coverage:

- Constructor restrictions (unsupported 2D theory combinations).
- `init`, `computeLfc`, and `initialGuessFromInput` behavior.
- Utility `Slfc` 2D/3D integration branches, `x=0` special case.
- Bridge-function and interpolation coupling behavior.

## 7) `esa.hpp` + `ESAUtil`

Classes:

- `ESA`
- `ESAUtil::Slfc`

Coverage:

- ESA top-level `compute` and resulting LFC shape/value sanity.
- `Slfc::get`, coefficient caching (`valid` transitions), and repeated-call consistency.
- NN/CSR/activation/free-energy paths with finite-value checks.
- Ground-state branch behavior in free-energy helper.

## 8) `qstls.hpp` + `QstlsUtil`

Classes:

- `Qstls`
- quantum utility classes in `QstlsUtil` (ADR/SSF and helpers)

Coverage:

- Constructor guards (unsupported dimensions/state points).
- Minimal successful compute path with tiny grids.
- ADR/SSF utility object `get()` paths including edge values.
- Persistence helpers (fixed components load/save) and error branches.

## 9) `qstlsiet.hpp` + `QstlsIetUtil`

Classes:

- `QstlsIet`
- utility SLFC/bridge integration helpers

Coverage:

- Constructor/compute constraints and mapping-dependent behavior.
- Utility integration branch behavior and finite-output assertions.
- Interaction with fixed data / guess paths.

## 10) `vsstls.hpp` + `qvsstls.hpp`

Classes:

- `VSStlsWorker`, `VSStlsManager`, `VSStls`
- `VSQstlsWorker`, `VSQstlsManager`, `QVSStls`
- `QAdder`

Coverage:

- Worker-level API contract (all getters + update/apply methods).
- Manager-level 3x3 orchestration behavior and derivative propagation.
- Scheme-level constructor guards and successful compute path (small grids).
- QAdder deterministic invariants (e.g., `S(q)=1` baseline).

## Execution order

1. `input` and deterministic util-heavy classes (`HFUtil`, `RpaUtil`).
2. STLS + IET families (`Stls`, `StlsIet`, `IetUtil`).
3. ESA branch and coefficient/caching tests.
4. Quantum families (`Qstls`, `QstlsIet`, QVS) with temp db/blob fixtures.
5. Final sweep for missing public methods and uncovered error paths.

## Quality gates for this module

- Every public method in `include/schemes/*.hpp` has at least one direct test.
- Every documented error branch has at least one failure-path test.
- Util namespace classes each have dedicated tests (not only indirect coverage).
- `./devtool build --native-only --native-tests` and `./devtool test native-cpp` pass.

## Execution environment note

- Run the validation commands inside the devcontainer, not on the host machine.
- Current container: `4deb4a17d562`.
- Example:
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool build --native-only --native-tests'`
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool test native-cpp'`
