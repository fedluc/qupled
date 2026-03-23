# Refactoring Plan – Step 3

## Goals

1. Rename `IVSWorker` → `VSWorker` (more conventional, matches project style).
2. Rename `StatePointGridBase` → `VSManager` (clarifies its role as coordinator).
3. Rename files `state_point_grid.{hpp,cpp}` → `vs_master.{hpp,cpp}`.
4. Move iteration method bodies from `VSStlsManager` / `VSQstlsManager`
   into `VSManager`, eliminating the duplicated `.cpp` boilerplate.

---

## Detailed changes

### 1. Rename `IVSWorker` → `VSWorker`

Replace all occurrences across:
- `include/vs/vsworker_base.hpp` (class definition + usages)
- `include/vs/vsstls.hpp` (`VSStlsWorker : public VSWorker, ...`)
- `include/vs/qvsstls.hpp` (`VSQstlsWorker : public VSWorker, ...`)

### 2. Rename `StatePointGridBase` → `VSManager`

Replace all occurrences across:
- `include/vs/vsmaster_base.hpp` (class definition)
- `include/vs/vsstls.hpp` (`VSStlsManager : public VSManager, ...`)
- `include/vs/qvsstls.hpp` (`VSQstlsManager : public VSManager, ...`)
- `src/vs/vsmaster_base.cpp` (all method definitions)
- `src/vs/vsstls.cpp` (constructor initialiser list, qualified calls)
- `src/vs/qvsstls.cpp` (constructor initialiser list, qualified calls)

### 3. Rename files

- `include/vs/vsworker_base.hpp` → `include/vs/vsworker.hpp`
- `include/vs/vsmaster_base.hpp` → `include/vs/vsmanager.hpp`
- `src/vs/vsmaster_base.cpp`     → `src/vs/vsmanager.cpp`

Update `#include` directives in every file that includes these headers:
- `include/vs/vsstls.hpp`
- `include/vs/qvsstls.hpp`
- `src/vs/vsbase.cpp` (if included transitively via vsbase.hpp, check)

Update `src/CMakeLists.txt`: replace `vs/vsmaster_base.cpp` with `vs/vsmanager.cpp`.

### 4. Move iteration method bodies into `VSManager`

The following methods have identical bodies in both `VSStlsManager` and
`VSQstlsManager`. Because `VSManager` is not in the `Stls`/`HF` inheritance
hierarchy, C++ cannot automatically use its methods to satisfy `HF`'s pure virtuals.
The overrides in the subclasses therefore cannot be removed, but their bodies can be
moved to `VSManager` as non-virtual protected helpers and the subclass overrides
reduced to inline one-liners in the headers.

Methods to move into `VSManager` (as non-virtual protected):

```cpp
// include/vs/vsmanager.hpp  (inside VSManager, protected section)
void init() {
  if (initDone) return;
  for (auto &w : workers) w->init();
  initDone = true;
}
void computeLfc()       { computeSynchronizedLfc(); }
void computeSsf()       { for (auto &w : workers) w->computeSsf(); }
double computeError() const {
  lastError = workers[GridPoints::CENTER.toIndex()]->computeError();
  return lastError;
}
void updateSolution()   { for (auto &w : workers) w->updateSolution(); }
void initialGuess()     { for (auto &w : workers) w->initialGuess(); }
```

Subclass headers become inline one-liners (no bodies in `.cpp`):

```cpp
// VSStlsManager / VSQstlsManager
void   init()                  override { VSManager::init(); }
void   computeLfc()            override { VSManager::computeLfc(); }
void   computeSsf()            override { VSManager::computeSsf(); }
double computeError()    const override { return VSManager::computeError(); }
void   updateSolution()        override { VSManager::updateSolution(); }
void   initialGuess()          override { VSManager::initialGuess(); }
```

Remove the corresponding method bodies from `src/vs/vsstls.cpp` and
`src/vs/qvsstls.cpp`.

---

## Files affected

| File | Change |
|------|--------|
| `include/vs/vsworker_base.hpp` | Rename → `vsworker.hpp`; rename class to `VSWorker` |
| `include/vs/vsmaster_base.hpp` | Rename → `vsmanager.hpp`; rename class to `VSManager` |
| `include/vs/vsstls.hpp` | Update include, rename bases to `VSManager`, rename master to `VSStlsManager`, inline overrides |
| `include/vs/qvsstls.hpp` | Update include, rename bases to `VSManager`, rename master to `VSQstlsManager`, inline overrides |
| `src/vs/vsmaster_base.cpp` | Rename → `vsmanager.cpp`; rename class to `VSManager` |
| `src/vs/vsstls.cpp` | Rename constructor to `VSStlsManager`, update qualified calls |
| `src/vs/qvsstls.cpp` | Rename constructor to `VSQstlsManager`, update qualified calls |
| `src/CMakeLists.txt` | Update source file name |

---

## What is NOT changed in this step

- `QAdder` class and its private default constructor (kept as-is).
- `computeQRaw` virtual method (kept; it is the Template Method hook in `VSBase`).
- `applyLfcDiff` / `computeBaseLfc` in worker classes (kept; they are pure virtual
  overrides from `VSWorker`).
- Public API of `VSStls` and `QVSStls` (unchanged).
