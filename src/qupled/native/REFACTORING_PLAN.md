# VS Code Refactoring Plan

## Goals

Replace the current CSR manager/worker pattern and `ThermoPropBase` hierarchy with a
cleaner composition-based design that separates three distinct responsibilities:

1. **State point coordination** (`StatePointGrid`) — creates and coordinates the 9
   scheme instances, handles the synchronized LFC calculation
2. **Thermodynamic integration and alpha solving** (`VSBase`) — absorbs all
   `ThermoPropBase` logic, owns the unified `computeAlpha()` formula
3. **Scheme orchestration** (`VSStls` / `QVSStls`) — single inheritance from `VSBase`,
   use composition with `StatePointGrid`

---

## Files Changed

### New files
- `state_point_grid.hpp`
- `state_point_grid.cpp`

### Significantly refactored
- `vsbase.hpp` / `vsbase.cpp`
- `vsstls.hpp` / `vsstls.cpp`
- `qvsstls.hpp` / `qvsstls.cpp`

### No existing files modified
`HF`, `Stls`, and `Qstls` are left completely untouched. All VS-specific behaviour is
contained in the new worker wrapper classes (see below).

---

## Classes Eliminated

| Class | Replacement |
|-------|-------------|
| `ThermoPropBase` | Logic absorbed into `VSBase` |
| `ThermoProp` | Eliminated (was only a constructor) |
| `QThermoProp` | Eliminated (logic moves to `QVSStls::computeQData()`) |
| `CSR` | Replaced by `StatePointGrid<Scheme>` |
| `StlsCSR` | Eliminated (replaced by `VSStlsWorker`) |
| `QstlsCSR` | Eliminated (replaced by `VSQstlsWorker`) |
| `StructIdx` enum | Replaced by `GridPoint` strong type |
| `ThermoIdx` enum | Internalized as fixed array indices in `VSBase` |

---

## New / Redesigned Classes

---

### `GridPoint` (in `state_point_grid.hpp`)

Strong type replacing the flat `StructIdx` / `ThermoIdx` enums.

```cpp
struct GridPoint {
    enum class Rs    { DOWN = -1, CENTER = 0, UP = 1 };
    enum class Theta { DOWN = -1, CENTER = 0, UP = 1 };

    Rs    rs;
    Theta theta;

    // Maps to a flat 0–8 index (theta-outer, rs-inner order)
    size_t toIndex() const;

    // Named constants for all 9 points
    static const GridPoint RS_DOWN_THETA_DOWN;
    static const GridPoint RS_THETA_DOWN;
    static const GridPoint RS_UP_THETA_DOWN;
    static const GridPoint RS_DOWN_THETA;
    static const GridPoint CENTER;            // {Rs::CENTER, Theta::CENTER}
    static const GridPoint RS_UP_THETA;
    static const GridPoint RS_DOWN_THETA_UP;
    static const GridPoint RS_THETA_UP;
    static const GridPoint RS_UP_THETA_UP;
};
```

---

### `VSStlsWorker` (in `vsstls.hpp`)
### `VSQstlsWorker` (in `qvsstls.hpp`)

Thin wrapper classes that inherit from `Stls` / `Qstls` and add only the
VS-specific behaviour needed by `StatePointGrid`. This keeps all VS concerns out
of the general-purpose `HF`, `Stls`, and `Qstls` classes.

Because they inherit from `HF` (via the chain `HF ← Rpa ← Stls ← Qstls`), they have
direct protected access to `lfc`, `mu`, `wvg`, `ssf`, and `adrFixedDatabaseName` —
no getters or setters need to be added to any existing class.

```cpp
// in vsstls.hpp
class VSStlsWorker : public Stls {
public:
    using Stls::Stls; // inherit constructors directly
    void applyLfcDiff(const Vector2D& v) { lfc.diff(v); }
};

// in qvsstls.hpp
class VSQstlsWorker : public Qstls {
public:
    // dbName is set here rather than needing a setter on Qstls
    VSQstlsWorker(const std::shared_ptr<const QVSStlsInput>& in,
                  const std::string& dbName);
    void   applyLfcDiff(const Vector2D& v) { lfc.diff(v); }
    // Builds QAdder using internal mu, wvg, ssf — mu never needs a public getter
    double computeQAdder(const std::shared_ptr<Integrator2D>& itg2D,
                         const std::vector<double>& itgGrid) const;
};
```

`StatePointGrid` is then instantiated with these worker types:
```cpp
// state_point_grid.cpp
template class StatePointGrid<VSStlsWorker>;
template class StatePointGrid<VSQstlsWorker>;
```

---

### `StatePointGrid<Scheme>` (in `state_point_grid.hpp` / `.cpp`)

Template class owning the 9 worker instances and coordinating their LFC calculation.
The manager/worker ambiguity of `CSR` is gone — this class is *only* the manager.

**Internal derivative tracking** — replaces the raw-pointer `DerivativeData` in `CSR`.
Indices into the workers array are used instead of raw pointers to neighbor LFCs:

```cpp
struct DerivativeData {
    enum class Type { CENTERED, FORWARD, BACKWARD };
    Type   type;
    size_t upIdx;
    size_t downIdx;
};
```

**The 3-step synchronized LFC** (preserves the critical ordering from `CSR::computeLfc()`):

```
Step 1: for each worker → w.computeLfc()             (base LFC, all workers first)
Step 2: for each worker → computeLfcDerivative(i)    (uses neighbor LFCs, safe now)
Step 3: for each worker → w.applyLfcDiff(lfcDerivatives[i])
```

`applyLfcDiff` is defined on `VSStlsWorker` / `VSQstlsWorker`, not on `HF` or `Stls`.

**Public interface:**

```cpp
template<typename Scheme>
class StatePointGrid {
public:
    StatePointGrid(/* scheme-specific input */, const VSInput& vsIn);

    void setAlpha(double alpha);
    int  compute();

    // Getters by GridPoint
    const std::vector<double>& getSsf(GridPoint p)  const;
    const Vector2D&            getLfc(GridPoint p)  const;
    const std::vector<double>& getWvg(GridPoint p)  const;
    double getCoupling(GridPoint p)          const;
    double getDegeneracy(GridPoint p)        const;
    double getUInt(GridPoint p)              const;
    double getFxcIntegrandValue(GridPoint p) const;
    double getAlpha()  const;   // central point
    double getError()  const;   // central point

    // Direct access to central worker (for idr, sdr, wvg delegation)
    const Scheme& centralWorker() const;

private:
    static constexpr int N = 9;
    std::array<std::unique_ptr<Scheme>, N> workers;
    std::array<Vector2D, N>                lfcDerivatives;
    double                                 alpha;

    std::array<DerivativeData, N> rsDerivData;
    std::array<DerivativeData, N> thetaDerivData;

    void   setupDerivativeData(const VSInput& vsIn);
    void   computeLfc();
    void   computeLfcDerivatives();
    double derivative(const Vector2D& f, int l, size_t i, DerivativeData::Type t) const;
    double derivative(double f0, double f1, double f2, DerivativeData::Type t)    const;
};
```

**Note on `StatePointGrid<VSQstlsWorker>` initialization** — the correct
`adrFixedDatabaseName` is passed to each `VSQstlsWorker` at construction time (not
via a setter), based on the worker's theta offset:

| Theta offset | Database name |
|---|---|
| `DOWN` | `{theory}_THETA_DOWN` |
| `CENTER` | `{theory}_THETA` |
| `UP` | `{theory}_THETA_UP` |

---

### `VSBase` (simplified, absorbs `ThermoPropBase`)

Single responsibility: orchestrate the secant solver and free energy integration.
The `thermoProp` shared pointer and all `ThermoPropBase` forwarding are removed.

**Key change: `computeAlpha()` is no longer virtual.**
The two alpha formulas are unified since the classical Q-term reduces to the internal
energy. `computeQData()` is the only virtual hook needed:

```
numerator   = Q + (1/3)·fxcr − (1/6)·fxcrr − (2/3)·(fxctt + fxcrt) + (1/3)·fxct
denominator = Q + (1/3)·Qr  + (2/3)·Qt
```

where for `VSStls`:  `{Q, Qr, Qt}` = `{uint, uintr, uintt}`
and for `QVSStls`: `{Q, Qr, Qt}` = quantum Q-term and its derivatives.

**Interface:**

```cpp
class VSBase : public Logger {
public:
    int    compute();
    double getAlpha() const { return alpha; }
    const std::vector<std::vector<double>>& getFreeEnergyIntegrand() const;
    const std::vector<double>&              getFreeEnergyGrid()      const;

protected:
    double alpha;
    // From ThermoPropBase (now owned here)
    std::vector<double>              rsGrid;
    std::vector<std::vector<double>> fxcIntegrand; // [3][nrs]: DOWN/CENTER/UP in theta

    // Abstract interface — must be implemented by VSStls / QVSStls
    virtual const VSInput& in() const = 0;
    virtual void init()            = 0;
    virtual void updateSolution()  = 0;

    // Grid access abstraction (delegates to StatePointGrid in derived classes)
    virtual int    runGrid()                                = 0;
    virtual double getCoupling(GridPoint p)    const        = 0;
    virtual double getDegeneracy(GridPoint p)  const        = 0;
    virtual double getFxcIntegrandValue(GridPoint p) const  = 0;

    // {Q, Qr, Qt}: classical = {uint, uintr, uintt}; quantum = Q-term derivatives
    virtual std::vector<double> computeQData() = 0;

    // Non-virtual helpers (from ThermoPropBase)
    void   setRsGrid();
    void   setFxcIntegrand();           // loads from input if provided
    void   updateFxcIntegrand();        // called after each runGrid()
    double computeFreeEnergy(GridPoint p, bool normalize) const;
    std::vector<double> getFreeEnergyData() const;  // {fxc, fxcr, fxcrr, fxct, fxctt, fxcrt}
    GridPoint getOutputGridPoint() const;           // handles rs=0 / theta=0 edge cases

private:
    // Unified — not virtual
    double computeAlpha();
    void   doIterations();
    double alphaDifference(const double& alphaTmp);
};
```

**`getOutputGridPoint()`** preserves the current `ThermoPropBase::getStructPropIdx()` logic:
- rs=0, theta=0  → `GridPoint::RS_DOWN_THETA_DOWN`
- rs>0, theta=0  → `GridPoint::RS_THETA_DOWN`
- rs=0, theta>0  → `GridPoint::RS_DOWN_THETA`
- rs>0, theta>0  → `GridPoint::CENTER`

---

### `VSStls` (single inheritance from `VSBase`, composition)

No longer inherits from `Stls`. Delegates all structural property access to the
internal `StatePointGrid<VSStlsWorker>`.

```cpp
class VSStls : public VSBase {
public:
    explicit VSStls(const std::shared_ptr<const VSStlsInput>& in);
    using VSBase::compute;

    // Public output interface (required by Python bindings)
    const std::vector<double>& getSsf()  const { return ssf; }
    const Vector2D&            getLfc()  const { return lfc; }
    const std::vector<double>& getWvg()  const;
    const Vector2D&            getIdr()  const;  // delegates to grid.centralWorker()
    std::vector<double>        getSdr()  const;  // delegates to grid.centralWorker()
    double                     getUInt() const;
    double                     getError() const; // delegates to grid.getError()

private:
    std::shared_ptr<const VSStlsInput> inPtr;
    StatePointGrid<VSStlsWorker>       grid;
    // Output storage (central state point, populated in updateSolution)
    std::vector<double> ssf;
    Vector2D            lfc;

    const VSInput& in() const override;
    void init()           override;
    void updateSolution() override; // copies ssf/lfc from grid at getOutputGridPoint()

    // VSBase abstract implementations
    int    runGrid()                               override;
    double getCoupling(GridPoint p)    const override;
    double getDegeneracy(GridPoint p)  const override;
    double getFxcIntegrandValue(GridPoint p) const override;
    // Returns {uint, uintr, uintt} — current ThermoPropBase::getInternalEnergyData() logic
    std::vector<double> computeQData() override;
};
```

---

### `QVSStls` (single inheritance from `VSBase`, composition)

Same pattern as `VSStls`. Additionally owns `itg2D` for QAdder computation.
`QAdder` class itself is unchanged.

```cpp
class QVSStls : public VSBase {
public:
    explicit QVSStls(const std::shared_ptr<const QVSStlsInput>& in);
    using VSBase::compute;

    // Same public interface as VSStls
    const std::vector<double>& getSsf()  const { return ssf; }
    const Vector2D&            getLfc()  const { return lfc; }
    const std::vector<double>& getWvg()  const;
    const Vector2D&            getIdr()  const;
    std::vector<double>        getSdr()  const;
    double                     getUInt() const;
    double                     getError() const;

private:
    std::shared_ptr<const QVSStlsInput> inPtr;
    StatePointGrid<VSQstlsWorker>       grid;
    std::shared_ptr<Integrator2D>       itg2D;   // for QAdder computation
    std::vector<double>                 itgGrid; // for segregated integration
    std::vector<double> ssf;
    Vector2D            lfc;

    const VSInput& in() const override;
    void init()           override;
    void updateSolution() override;

    int    runGrid()                               override;
    double getCoupling(GridPoint p)    const override;
    double getDegeneracy(GridPoint p)  const override;
    double getFxcIntegrandValue(GridPoint p) const override;
    // Returns {Q, Qr, Qt} using VSQstlsWorker::computeQAdder() on relevant grid workers
    // Current QThermoProp::getQData() logic, moved here
    std::vector<double> computeQData() override;
};
```

---

## Python Bindings

`schemes.cpp` requires no changes. `PyScheme<VSStls, VSStlsInput>` and
`PyScheme<QVSStls, QVSStlsInput>` will continue to work because `VSStls` and `QVSStls`
provide all methods called by `exposeVSSchemeClass` and `exposeIterativeSchemeProperties`:

| Python property | Method called | Source in new design |
|---|---|---|
| `compute` | `VSStls::compute()` | inherited from `VSBase` |
| `lfc` | `VSStls::getLfc()` | output storage, set in `updateSolution()` |
| `ssf` | `VSStls::getSsf()` | output storage |
| `wvg` | `VSStls::getWvg()` | delegates to `grid.centralWorker()` |
| `idr` | `VSStls::getIdr()` | delegates to `grid.centralWorker()` |
| `sdr` | `VSStls::getSdr()` | delegates to `grid.centralWorker()` |
| `uint` | `VSStls::getUInt()` | computed from stored ssf/wvg/rs |
| `error` | `VSStls::getError()` | delegates to `grid.getError()` |
| `alpha` | `VSBase::getAlpha()` | inherited from `VSBase` |
| `free_energy_integrand` | `VSBase::getFreeEnergyIntegrand()` | inherited from `VSBase` |
| `free_energy_grid` | `VSBase::getFreeEnergyGrid()` | inherited from `VSBase` |

---

## Responsibility Summary

| Class | Responsibility |
|-------|---------------|
| `VSBase` | Secant solver loop, rs grid, fxc integrand storage, free energy data, unified alpha formula |
| `StatePointGrid<T>` | 9-point construction and coordination, synchronized 3-step LFC, getters by `GridPoint` |
| `VSStlsWorker` | `Stls` + `applyLfcDiff` — VS-specific worker for `StatePointGrid<VSStlsWorker>` |
| `VSQstlsWorker` | `Qstls` + `applyLfcDiff` + `computeQAdder` — VS-specific worker for `StatePointGrid<VSQstlsWorker>` |
| `VSStls` | Wires `VSBase` + `StatePointGrid<VSStlsWorker>`, computes `{uint, uintr, uintt}`, exposes output |
| `QVSStls` | Wires `VSBase` + `StatePointGrid<VSQstlsWorker>`, computes quantum `{Q, Qr, Qt}`, exposes output |
| `HF` / `Stls` / `Qstls` | Completely unchanged |
| `GridPoint` | Strong-typed address into the 3×3 state point grid |
| `QAdder` | Unchanged |
