# C++ Testing Guidelines (Reusable Blueprint)

## Purpose

Use this document as a cross-module blueprint for writing and reviewing native C++ tests.
It is distilled from the Schemes and Thermo test plans.

## 1) Core principles

- Test public behavior, not implementation details.
- Cover every public method in scope at least once.
- Treat utility helpers/classes as first-class test targets, not only indirectly covered code.
- Keep tests deterministic: small grids, explicit tolerances, no hidden randomness.
- Assert both value correctness and contract/invariant correctness (shape, dimensions, sizes, sentinel states).
- Test both success paths and failure paths.

## 2) Granularity rules (important)

- Default to one test per behavior/branch.
- Prefer one test for each branch split:
  - 2D vs 3D
  - `tau = 0` vs `tau > 0`
  - `x = 0` vs `x > 0`
  - zero coupling vs finite coupling
- Keep one clear failure reason per test.
- If a single test constructs multiple independent objects/branches, treat that as a refactor signal and split it.
- Group assertions only when they represent one inseparable behavior and share expensive setup.

## 3) Constructor coverage

- Add dedicated constructor tests when constructors:
  - enforce constraints,
  - initialize critical storage/state,
  - or branch on input configuration.
- For guarded constructors, include both:
  - `RejectsInvalid...`
  - `AcceptsValid...`
- Keep constructor tests separate from compute/algorithm tests.

## 4) Private methods policy

- Do not test private methods directly.
- Cover private logic through public observable behavior.
- If private logic is complex and hard to exercise, extract a helper and test that helper.
- Avoid visibility hacks (e.g. redefining `private`).

## 5) Test organization and naming

- Split test files by domain/class area for readability.
- Prefer descriptive names tied to behavior, e.g.:
  - `ReturnsZeroAtXZero`
  - `ValidatesLfcShapeForInteractingTheory`
  - `ComputesNormalizedIntegral`
- Keep test names specific enough that failures are self-explanatory.

## 6) Coverage planning workflow

For each header/module under test:

1. List public classes/functions in scope.
2. Enumerate behavior branches (dimension/state/mode/error cases).
3. Map each branch to at least one dedicated test.
4. Add edge/boundary checks and known analytic/reference points when available.
5. Add explicit error-path assertions (`EXPECT_THROW`, invalid shape/value guards).

## 7) Operational checklist (easy to overlook)

## Validation environment

- Run validation inside the devcontainer, not the host machine.
- Current container used in this repo:
  - `4deb4a17d562`
- Validation commands:
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool build --native-only --native-tests'`
  - `docker exec 4deb4a17d562 /bin/bash -lc 'cd /workspaces/qupled && ./devtool test native-cpp'`

## Commit author reminder

- For test-plan/test-execution commits in this workflow, set author to:
  - `Federico Lucco Castello <71817492+fedluc@users.noreply.github.com>`

Example:

- `git commit --author="Federico Lucco Castello <71817492+fedluc@users.noreply.github.com>" -m "<message>"`

## 8) Definition of done

- Public methods/functions in scope are covered.
- Branch/error coverage is explicit and modular.
- Tests remain readable and diagnosable.
- Native build and test validation pass in devcontainer.
- Commit metadata (author) is correct when required.
