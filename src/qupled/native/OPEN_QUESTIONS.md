# Open Questions for Step 2 Refactoring

## Background

Two improvements are being planned for `src/qupled/native`:

1. Replace the template-based `StatePointGrid<Scheme>` with dynamic dispatch /
   virtual functions, so that the iteration loop is not reimplemented and new
   schemes (e.g. VSRpa) can be added without forcing them into the STLS
   iteration structure.

2. Unify `QAdder` into a single class that handles both the classical case
   (Q = internal energy) and the quantum case, selected by a mode parameter.

---

## Issue 1: StatePointGrid — replacing templates with dynamic dispatch

### Proposed approach (to confirm)

1. Define an abstract `IVSWorker` interface (base class for all workers), exposing
   only what `StatePointGrid` needs to orchestrate the synchronized LFC:
   `computeBaseLfc()`, `applyLfcDiff()`, `getLfc()`, `getWvg()`, and value getters.
2. Make `StatePointGrid` a non-template base class using
   `std::array<std::unique_ptr<IVSWorker>, N>`.
3. Make `StatePointGrid::compute()` pure virtual so that concrete subclasses
   can provide different iteration structures:
   - `StatePointGridIterative` — for STLS/Qstls-style iterative schemes
   - `StatePointGridDirect` (future) — for non-iterative schemes like RPA
4. `StatePointGridIterative` defines a richer `IIterativeVSWorker : IVSWorker`
   sub-interface that adds step methods (`doInitialGuess`, `doComputeSsf`, etc.).

### Questions

- **Goal clarity:** Is the main goal to make adding new scheme types (e.g. VSRpa)
  easier, or to literally eliminate the reimplemented loop? If VSRpa is
  non-iterative, `StatePointGridDirect::compute()` would have a completely
  different body anyway — so the iteration is not eliminated, just moved to a
  subclass.

- **VSRpa iteration structure:** What does VSRpa's iteration look like exactly?
  Does it have any convergence loop at all, or is it a single-shot computation?
  This determines whether a `StatePointGridDirect` makes sense or something more
  specific is needed.

- **Worker step granularity:** The current worker `do*()` methods are thin wrappers
  around protected `Stls` methods (`doComputeSsf` → `Stls::computeSsf()` etc.).
  In the virtual dispatch model these become virtual methods on
  `IIterativeVSWorker`. Should we keep this level of granularity (separate
  virtual calls for init / initialGuess / computeSsf / computeError /
  updateSolution), or collapse them into fewer calls like `doOneIterationStep()`?

---

## Issue 2: Unified QAdder

### Proposed approach (to confirm)

- Single `QAdder` class with a mode parameter (`CLASSICAL` / `QUANTUM`).
- `CLASSICAL` mode: computes the internal energy (currently done directly in
  `VSStls::computeQData()` via `thermoUtil::computeInternalEnergy`).
- `QUANTUM` mode: runs the current 2D-integral formula (current `QAdder::get()`).
- Both `VSStls::computeQData()` and `QVSStls::computeQData()` use this unified class.

### Questions

- **Classical mode implementation:** Should `QAdder::get()` in classical mode
  delegate to `thermoUtil::computeInternalEnergy` (simple but heterogeneous), or
  recompute via the 1D integrator itself (uniform but redundant)? The former is
  simpler; the latter makes the two modes more symmetric.

- **Location:** The class currently lives in `qvsstls.hpp`. Since both `VSStls`
  and `QVSStls` would use it after unification, it probably belongs in
  `vsbase.hpp` or its own small header. What is the preference?

- **Can `computeQData()` become non-virtual?** It is currently pure virtual in
  `VSBase` — the only difference between `VSStls` and `QVSStls` is which mode
  is used. If `QAdder` is unified, could `computeQData()` itself become
  non-virtual in `VSBase`, with the mode selection being the only
  subclass-specific part? Or do the two implementations differ enough in how
  they access the grid that keeping it virtual makes more sense?
