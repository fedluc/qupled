# VS Code Refactoring — Step 2

## Goals

1. **Replace the `StatePointGrid<Scheme>` template with a class hierarchy** so that
   - the iteration loop is no longer re-implemented inside `StatePointGrid`,
   - new VS schemes (e.g. `VSRpa`) can be added without being forced into the STLS
     step-by-step iteration structure.

2. **Unify `QAdder`** into a single class (classical and quantum modes) placed in
   `vsbase.hpp`, and make `VSBase::computeQData()` non-virtual.

---

## Issue 1 — Class hierarchy instead of template

### Key insight

`HF::compute()` already calls the virtual `computeStructuralProperties()`.
`Stls` overrides that method with the full iteration loop, which in turn calls
the virtual chain `computeLfc()` → `computeSsf()` → `computeError()` →
`updateSolution()`.

A class `StatePointGridVSStls` that inherits from both `StatePointGridBase`
and `Stls` only needs to override those individual virtual methods to coordinate
all 9 workers. **The loop itself is inherited from `Stls` at no cost** — no
re-implementation required.

The same pattern applies to `StatePointGridVSQstls` (inherits from `Qstls`) and
to any future `StatePointGridVSRpa` (would inherit from the appropriate base).

### New / modified classes

---

#### `IVSWorker`  (new, in `state_point_grid.hpp`)

Abstract interface for all worker objects. Exposes:

```
// For StatePointGridBase — synchronized LFC
computeBaseLfc()
applyLfcDiff(const Vector2D&)

// Getters used by StatePointGridBase
getLfc(), getWvg(), getSsf()
getCoupling(), getDegeneracy()
getUInt(), getFxcIntegrandValue()

// For StatePointGridVSStls/QVSStls iteration coordination
workerInit()
workerInitialGuess()
workerComputeSsf()
workerComputeError() → double
workerUpdateSolution()
```

The `worker*` methods allow `StatePointGridVSStls` to drive each worker's
individual iteration steps without needing to know the concrete type.
Non-iterative future workers (e.g. `VSRpaWorker`) would leave `workerInitialGuess`,
`workerComputeError`, etc. as no-ops.

---

#### `StatePointGridBase`  (replaces non-template `StatePointGrid`, in `state_point_grid.hpp/cpp`)

Manages the 9 worker instances and owns all derivative/grid bookkeeping.
Identical responsibilities to the current `StatePointGrid`, minus the iteration loop.

```cpp
class StatePointGridBase {
public:
    void setAlpha(double);
    double getAlpha() const;
    double getError() const;

    // Getters by GridPoint — delegate to workers array
    const vector<double>& getSsf(GridPoint) const;
    const Vector2D&        getLfc(GridPoint) const;
    const vector<double>& getWvg(GridPoint) const;
    double getCoupling(GridPoint)         const;
    double getDegeneracy(GridPoint)       const;
    double getUInt(GridPoint)             const;
    double getFxcIntegrandValue(GridPoint)const;
    const IVSWorker& getWorkerAt(GridPoint) const;

protected:
    static constexpr int N = 9;
    array<unique_ptr<IVSWorker>, N> workers;
    array<Vector2D, N>              lfcDerivatives;
    array<DerivativeData, N>        rsDerivData;
    array<DerivativeData, N>        thetaDerivData;
    array<double, N>                rsValues;
    array<double, N>                thetaValues;
    double alpha;
    double lastError;

    // Called by StatePointGridVSStls::computeLfc()
    void computeSynchronizedLfc();

private:
    void   setupDerivativeData();
    void   computeLfcDerivatives();
    double derivative(...) const;  // two overloads, unchanged
};
```

The constructor is protected/delegated — workers are populated by the derived
class before `setupDerivativeData()` is called.

---

#### `StatePointGridVSStls`  (new, in `vsstls.hpp/cpp`)

```cpp
class StatePointGridVSStls : public StatePointGridBase, public Stls {
public:
    explicit StatePointGridVSStls(const shared_ptr<const VSStlsInput>& in,
                                  bool verbose = false);
    // compute() is inherited from HF (calls init() then computeStructuralProperties())

protected:
    // Override Stls/HF virtual methods — each fans out to all 9 workers
    void   init()             override;  // workers[i]->workerInit() for all i
    void   computeLfc()       override;  // calls computeSynchronizedLfc()
    void   computeSsf()       override;  // workers[i]->workerComputeSsf()
    double computeError()const override; // workers[CENTER]->workerComputeError()
    void   updateSolution()   override;  // workers[i]->workerUpdateSolution()
    void   initialGuess()     override;  // workers[i]->workerInitialGuess()
};
```

Constructor:
1. Initialise `Stls` with the central input and `verbose=false` (suppresses
   per-iteration logging; the outer VSStls prints the summary itself).
2. Build 9 `VSStlsWorker` instances from perturbed inputs (same arithmetic as
   the current `StatePointGrid` template constructor).
3. Populate `workers[]`, `rsValues[]`, `thetaValues[]` in `StatePointGridBase`.
4. Call `setupDerivativeData()`.

**The 9 workers array includes the central point.** The `Stls` parent of
`StatePointGridVSStls` provides the iteration framework; the actual structural
data lives entirely in the workers. The parent's own `ssf`/`lfc` fields are
not used.

---

#### `StatePointGridVSQstls`  (new, in `qvsstls.hpp/cpp`)

Same pattern, inherits from `StatePointGridBase` and `Qstls`.
The `Qstls::init()` and `Qstls::computeLfc()` overrides in the Qstls chain are
already virtual, so the same override scheme applies.

---

#### Worker changes

**`VSStlsWorker`** (in `vsstls.hpp`):
- Remove `doInit()`, `doInitialGuess()`, `doComputeLfc()`, `doComputeSsf()`,
  `doUpdateSolution()`, `doComputeError()`.
- Inherit `IVSWorker` and implement all interface methods (thin wrappers that
  call the protected parent methods, identical to the old `do*` methods but
  now under the interface contract).
- Keep `applyLfcDiff()`.
- `computeBaseLfc()` → calls `Stls::computeLfc()` explicitly (avoids accidentally
  calling `StatePointGridVSStls::computeLfc()` if ever called on the grid itself).

**`VSQstlsWorker`** (in `qvsstls.hpp`):
- Same changes.
- Keep `computeQAdder()`.

---

### Files changed

| File | Change |
|------|--------|
| `state_point_grid.hpp` | Replace template class with `IVSWorker` interface + `StatePointGridBase` |
| `state_point_grid.cpp` | Remove template; implement `StatePointGridBase`; remove explicit instantiations |
| `vsstls.hpp` | Add `StatePointGridVSStls`; update `VSStlsWorker` to implement `IVSWorker` |
| `vsstls.cpp` | Implement `StatePointGridVSStls`; update `VSStls` to use new grid type |
| `qvsstls.hpp` | Add `StatePointGridVSQstls`; update `VSQstlsWorker` to implement `IVSWorker` |
| `qvsstls.cpp` | Same |

`HF`, `Stls`, `Qstls` — **no changes required**. All the virtual hooks we need
(`init`, `computeLfc`, `computeSsf`, `computeError`, `updateSolution`,
`initialGuess`, `computeStructuralProperties`) are already virtual in the
existing hierarchy.

---

## Issue 2 — Unified `QAdder`

### Design

`QAdder` gains a `Mode` enum and two named constructors. The classical mode
delegates to `thermoUtil::computeInternalEnergy`; the quantum mode runs the
current 2D-integral formula.

```cpp
// in vsbase.hpp

class QAdder {
public:
    enum class Mode { CLASSICAL, QUANTUM };

    // Classical mode — computes internal energy
    static QAdder classical(const vector<double>& wvg,
                            const vector<double>& ssf,
                            const shared_ptr<const Input>& in);

    // Quantum mode — current formula (constructor arguments unchanged)
    static QAdder quantum(double Theta, double mu,
                          double limitMin, double limitMax,
                          const vector<double>& itgGrid,
                          shared_ptr<Integrator1D> itg1,
                          shared_ptr<Integrator2D> itg2,
                          shared_ptr<Interpolator1D> interp);

    double get() const;

private:
    Mode mode;
    // Classical fields
    vector<double> wvg_;
    vector<double> ssf_;
    shared_ptr<const Input> in_;
    // Quantum fields (current QAdder members)
    double Theta, mu;
    pair<double,double> limits;
    const vector<double>* itgGrid_;
    shared_ptr<Integrator1D>    itg1_;
    shared_ptr<Integrator2D>    itg2_;
    shared_ptr<Interpolator1D>  interp_;

    // Existing private methods (quantum path only)
    double ssf(const double& y) const;
    double integrandDenominator(double q) const;
    // ...
};
```

### Making `computeQData()` non-virtual

`VSBase` gets a new protected virtual hook `computeQRaw(GridPoint) → double`,
which returns Q *without* rs-normalisation (matching the existing asymmetry
between `getFxcIntegrandValue` — rs=1 — and the derivative formulas).

`VSBase::computeQData()` becomes a non-virtual concrete method that calls
`computeQRaw()` for the needed grid points and assembles the derivative vector:

```cpp
// non-virtual, in vsbase.cpp
vector<double> VSBase::computeQData() {
    const double rs  = getCoupling(CENTER);
    const double q   = computeQRaw(CENTER) / rs;
    const double drs = getCoupling(RS_UP_THETA) - rs;
    const double qr  = (computeQRaw(RS_UP_THETA) - computeQRaw(RS_DOWN_THETA))
                       / (2.0 * drs) - q;
    const double theta = getDegeneracy(CENTER);
    const double dt    = getDegeneracy(RS_THETA_UP) - theta;
    const double qt    = theta
                         * (computeQRaw(RS_THETA_UP)/rs - computeQRaw(RS_THETA_DOWN)/rs)
                         / (2.0 * dt);
    return {q, qr, qt};
}
```

`computeQRaw(GridPoint)` is the new (smaller) virtual surface:

```cpp
// VSStls — classical
double VSStls::computeQRaw(GridPoint p) const {
    return QAdder::classical(grid.getWvg(p), grid.getSsf(p), inPtr).get();
    // classical QAdder calls thermoUtil::computeInternalEnergy with rs=1
}

// QVSStls — quantum
double QVSStls::computeQRaw(GridPoint p) const {
    return grid.getWorkerAt(p).computeQAdder(itg2D, itgGrid);
}
```

### Files changed

| File | Change |
|------|--------|
| `vsbase.hpp` | Move `QAdder` here (from `qvsstls.hpp`); add `Mode`; add named constructors; add `virtual computeQRaw(GridPoint) = 0`; make `computeQData()` non-virtual |
| `vsbase.cpp` | Implement non-virtual `computeQData()`; implement `QAdder` classical/quantum paths |
| `vsstls.hpp` | Add `computeQRaw(GridPoint)` override; remove old `computeQData()` |
| `vsstls.cpp` | Implement `computeQRaw()` using `QAdder::classical` |
| `qvsstls.hpp` | Remove `QAdder` class definition; add `computeQRaw(GridPoint)` override; remove `computeQData()` |
| `qvsstls.cpp` | Implement `computeQRaw()` using `QAdder::quantum` |

---

---

## Issue 3 — Move VS code to a dedicated submodule

All VS-related headers and sources move into `include/vs/` and `src/vs/`
subdirectories to separate them from the non-VS code.

### New directory layout

```
include/
  vs/
    state_point_grid.hpp   ← was include/state_point_grid.hpp
    vsbase.hpp             ← was include/vsbase.hpp
    vsstls.hpp             ← was include/vsstls.hpp
    qvsstls.hpp            ← was include/qvsstls.hpp
src/
  vs/
    state_point_grid.cpp   ← was src/state_point_grid.cpp
    vsbase.cpp             ← was src/vsbase.cpp
    vsstls.cpp             ← was src/vsstls.cpp
    qvsstls.cpp            ← was src/qvsstls.cpp
```

### Include path convention

All `#include` directives that reference VS headers use the `vs/` prefix
from the common `include/` root (which stays as the only include directory):

```cpp
#include "vs/state_point_grid.hpp"  // not "state_point_grid.hpp"
#include "vs/vsbase.hpp"
// etc.
```

This applies both to files outside the VS module (e.g. `schemes.cpp`) and to
files inside it (e.g. `vs/vsstls.hpp` including `vs/vsbase.hpp`).

### Files changed

| File | Change |
|------|--------|
| `src/CMakeLists.txt` | Replace `state_point_grid.cpp`, `vsbase.cpp`, `vsstls.cpp`, `qvsstls.cpp` with `vs/` prefixed paths |
| `src/python_interface/schemes.cpp` | Update `#include "qvsstls.hpp"` → `"vs/qvsstls.hpp"` etc. |
| Old flat VS files | Deleted |

No changes to the `target_include_directories` CMake directive — the single
`../include` entry already covers `include/vs/` transitively.

---

## Responsibility summary after step 2

| Class | Responsibility |
|-------|---------------|
| `IVSWorker` | Abstract interface: synchronized LFC hooks + iteration step hooks |
| `StatePointGridBase` | Manages 9 workers; synchronized LFC; getters by GridPoint |
| `StatePointGridVSStls` | Inherits Stls iteration loop; overrides virtual methods to fan out to workers |
| `StatePointGridVSQstls` | Same pattern for Qstls |
| `VSStlsWorker` | Implements `IVSWorker`; thin wrappers around protected Stls methods |
| `VSQstlsWorker` | Same for Qstls; keeps `computeQAdder()` |
| `VSBase` | Secant solver, rs grid, fxc integrand, non-virtual `computeQData()`, `QAdder` |
| `VSStls` | Wires `VSBase` + `StatePointGridVSStls`; implements `computeQRaw()` classically |
| `QVSStls` | Wires `VSBase` + `StatePointGridVSQstls`; implements `computeQRaw()` quantumly |
| `HF` / `Stls` / `Qstls` | Completely unchanged |
| `QAdder` | Unified classical/quantum Q computation; lives in `vsbase.hpp` |
