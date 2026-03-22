# Refactoring Plan – Step 3

## Goals

1. Rename `IVSWorker` → `VSWorkerBase` (more conventional, matches project style).
2. Rename `StatePointGridBase` → `VSMasterBase` (clarifies its role as coordinator).
3. Rename files `state_point_grid.{hpp,cpp}` → `vs_master.{hpp,cpp}`.
4. Move iteration method bodies from `StatePointGridVSStls` / `StatePointGridVSQstls`
   into `VSMasterBase`, eliminating the duplicated `.cpp` boilerplate.

---

## Detailed changes

### 1. Rename `IVSWorker` → `VSWorkerBase`

Replace all occurrences across:
- `include/vs/state_point_grid.hpp` (class definition + usages)
- `include/vs/vsstls.hpp` (`VSStlsWorker : public IVSWorker, ...`)
- `include/vs/qvsstls.hpp` (`VSQstlsWorker : public IVSWorker, ...`)

### 2. Rename `StatePointGridBase` → `VSMasterBase`

Replace all occurrences across:
- `include/vs/state_point_grid.hpp` (class definition)
- `include/vs/vsstls.hpp` (`StatePointGridVSStls : public StatePointGridBase, ...`)
- `include/vs/qvsstls.hpp` (`StatePointGridVSQstls : public StatePointGridBase, ...`)
- `src/vs/state_point_grid.cpp` (all method definitions)
- `src/vs/vsstls.cpp` (constructor initialiser list, qualified calls)
- `src/vs/qvsstls.cpp` (constructor initialiser list, qualified calls)

### 3. Rename files

- `include/vs/state_point_grid.hpp` → `include/vs/vs_master.hpp`
- `src/vs/state_point_grid.cpp`     → `src/vs/vs_master.cpp`

Update `#include` directives in every file that includes `"vs/state_point_grid.hpp"`:
- `include/vs/vsstls.hpp`
- `include/vs/qvsstls.hpp`
- `src/vs/vsbase.cpp` (if included transitively via vsbase.hpp, check)

Update `src/CMakeLists.txt`: replace `vs/state_point_grid.cpp` with `vs/vs_master.cpp`.

### 4. Move iteration method bodies into `VSMasterBase`

The following methods have identical bodies in both `StatePointGridVSStls` and
`StatePointGridVSQstls`. Because `VSMasterBase` is not in the `Stls`/`HF` inheritance
hierarchy, C++ cannot automatically use its methods to satisfy `HF`'s pure virtuals.
The overrides in the subclasses therefore cannot be removed, but their bodies can be
moved to `VSMasterBase` as non-virtual protected helpers and the subclass overrides
reduced to inline one-liners in the headers.

Methods to move into `VSMasterBase` (as non-virtual protected):

```cpp
// include/vs/vs_master.hpp  (inside VSMasterBase, protected section)
void masterInit() {
  if (initDone) return;
  for (auto &w : workers) w->workerInit();
  initDone = true;
}
void masterComputeLfc()       { computeSynchronizedLfc(); }
void masterComputeSsf()       { for (auto &w : workers) w->workerComputeSsf(); }
double masterComputeError() const {
  lastError = workers[GridPoints::CENTER.toIndex()]->workerComputeError();
  return lastError;
}
void masterUpdateSolution()   { for (auto &w : workers) w->workerUpdateSolution(); }
void masterInitialGuess()     { for (auto &w : workers) w->workerInitialGuess(); }
```

Subclass headers become inline one-liners (no bodies in `.cpp`):

```cpp
// StatePointGridVSStls / StatePointGridVSQstls
void   init()                  override { masterInit(); }
void   computeLfc()            override { masterComputeLfc(); }
void   computeSsf()            override { masterComputeSsf(); }
double computeError()    const override { return masterComputeError(); }
void   updateSolution()        override { masterUpdateSolution(); }
void   initialGuess()          override { masterInitialGuess(); }
```

Remove the corresponding method bodies from `src/vs/vsstls.cpp` and
`src/vs/qvsstls.cpp`.

---

## Files affected

| File | Change |
|------|--------|
| `include/vs/state_point_grid.hpp` | Rename → `vs_master.hpp`; rename classes |
| `include/vs/vsstls.hpp` | Update include, rename bases, inline overrides |
| `include/vs/qvsstls.hpp` | Update include, rename bases, inline overrides |
| `src/vs/state_point_grid.cpp` | Rename → `vs_master.cpp`; rename class |
| `src/vs/vsstls.cpp` | Remove 6 method bodies, update qualified calls |
| `src/vs/qvsstls.cpp` | Remove 6 method bodies, update qualified calls |
| `src/CMakeLists.txt` | Update source file name |

---

## What is NOT changed in this step

- `QAdder` class and its private default constructor (kept as-is).
- `computeQRaw` virtual method (kept; it is the Template Method hook in `VSBase`).
- `applyLfcDiff` / `computeBaseLfc` in worker classes (kept; they are pure virtual
  overrides from `VSWorkerBase`).
- Public API of `VSStls` and `QVSStls` (unchanged).
