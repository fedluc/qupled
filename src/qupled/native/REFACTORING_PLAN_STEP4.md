# Refactoring Plan – Step 4

## Goals

1. **Remove debug output** from `VSQstlsWorker::init()`.
2. **Introduce virtual `grid()` method** in `VSBase` to provide unified access to the manager.
3. **Eliminate code repetition** by moving duplicated methods from `VSStls`/`QVSStls` to `VSBase`.
4. **Add convenience methods** to manager classes to support public API without adding methods to `VSWorker`.

---

## Key Insight

Currently, both `VSStls` and `QVSStls` have nearly identical implementations for many methods that simply delegate to their `grid` member. By introducing a virtual `grid()` accessor in `VSBase`, we can:

1. Move all duplicated delegation methods to `VSBase` as non-virtual implementations
2. Access the manager through the `grid()` virtual method
3. Add type-specific convenience methods to the manager classes (not the worker interface) to handle downcasting

This pattern is similar to how `in()` works: the child class stores a `std::shared_ptr<const VSStlsInput> inPtr`, and the base class accesses it through a virtual method that performs the upcast.

---

## Changes to Classes

### 1. `VSQstlsWorker` — Remove debug output

**File**: `include/vs/qvsstls.hpp`

Remove the debug print statement from `init()`:

```cpp
// Before:
void init() override {
  std::cerr << "Initializing VSQstlsWorker..." << std::endl;
  Qstls::init();
}

// After:
void init() override { Qstls::init(); }
```

---

### 2. `VSBase` — Add virtual `grid()` method

**File**: `include/vs/vsbase.hpp`

Add to the protected section:

```cpp
protected:
  // Virtual grid accessor (pure virtual)
  virtual VSManager& grid() = 0;
  virtual const VSManager& grid() const = 0;
```

Move these methods from derived classes to `VSBase` as **non-virtual** implementations:

```cpp
protected:
  // Abstract interface — REMOVE these lines (moved to private or kept as before)
  // virtual int    runGrid() = 0;                                // ← DELETE
  // virtual double getCoupling(GridPoint p)    const        = 0; // ← DELETE
  // virtual double getDegeneracy(GridPoint p)  const        = 0; // ← DELETE
  // virtual double getFxcIntegrandValue(GridPoint p) const  = 0; // ← DELETE

  // Still needed as virtual (scheme-specific logic)
  virtual const VSInput &in() const = 0;
  virtual const Input &inScheme() const = 0;
  virtual double computeQRaw(GridPoint p) const = 0;
```

Add these as **non-virtual protected** methods (delegate to `grid()`):

```cpp
protected:
  // Grid coordination (non-virtual, uses grid())
  int runGrid();
  double getCoupling(GridPoint p) const;
  double getDegeneracy(GridPoint p) const;
  double getFxcIntegrandValue(GridPoint p) const;
```

Add these as **non-virtual public** methods (delegate to `grid()`):

```cpp
public:
  // Public output interface (non-virtual, uses grid())
  const std::vector<double>& getSsf() const;
  const Vector2D& getLfc() const;
  const std::vector<double>& getWvg() const;
  const Vector2D& getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  double getError() const;
```

**File**: `src/vs/vsbase.cpp`

Implement the delegating methods:

```cpp
int VSBase::runGrid() {
  grid().setAlpha(alpha);
  int status = grid().compute();
  println(formatUtil::format("Alpha = {:.5e}, Residual error "
                             "(structural properties) = {:.5e}",
                             grid().getAlpha(),
                             grid().getError()));
  return status;
}

double VSBase::getCoupling(GridPoint p) const {
  return grid().getCoupling(p);
}

double VSBase::getDegeneracy(GridPoint p) const {
  return grid().getDegeneracy(p);
}

double VSBase::getFxcIntegrandValue(GridPoint p) const {
  return grid().getFxcIntegrandValue(p);
}

const std::vector<double>& VSBase::getSsf() const {
  return grid().getSsf();
}

const Vector2D& VSBase::getLfc() const {
  return grid().getLfc();
}

const std::vector<double>& VSBase::getWvg() const {
  return grid().getWvg();
}

const Vector2D& VSBase::getIdr() const {
  return grid().getIdr();
}

std::vector<double> VSBase::getSdr() const {
  return grid().getSdr();
}

double VSBase::getUInt() const {
  return grid().getUInt();
}

double VSBase::getError() const {
  return grid().getError();
}
```

---

### 3. `VSManager` — Add convenience methods for common worker access

**File**: `include/vs/vsmanager.hpp`

Add to public section:

```cpp
public:
  // Convenience getters for central worker (delegates to getWorkerAt(CENTER))
  const std::vector<double>& getSsf() const;
  const Vector2D& getLfc() const;
  const std::vector<double>& getWvg() const;
```

**File**: `src/vs/vsmanager.cpp`

Implement:

```cpp
const std::vector<double>& VSManager::getSsf() const {
  return getWorkerAt(GridPoints::CENTER).getSsf();
}

const Vector2D& VSManager::getLfc() const {
  return getWorkerAt(GridPoints::CENTER).getLfc();
}

const std::vector<double>& VSManager::getWvg() const {
  return getWorkerAt(GridPoints::CENTER).getWvg();
}
```

---

### 4. `VSStlsManager` — Add scheme-specific convenience methods

**File**: `include/vs/vsstls.hpp`

Add to `VSStlsManager` public section:

```cpp
public:
  // Scheme-specific getters (downcast to Stls)
  const Vector2D& getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
```

**File**: `src/vs/vsstls.cpp`

Implement in `VSStlsManager`:

```cpp
const Vector2D& VSStlsManager::getIdr() const {
  return dynamic_cast<const Stls&>(getWorkerAt(GridPoints::CENTER)).getIdr();
}

std::vector<double> VSStlsManager::getSdr() const {
  return dynamic_cast<const Stls&>(getWorkerAt(GridPoints::CENTER)).getSdr();
}

double VSStlsManager::getUInt() const {
  return dynamic_cast<const Stls&>(getWorkerAt(GridPoints::CENTER)).getUInt();
}
```

---

### 5. `VSStls` — Implement `grid()` and remove duplicated methods

**File**: `include/vs/vsstls.hpp`

Remove these method declarations from `VSStls` (now inherited from `VSBase`):

```cpp
// DELETE these lines:
// const std::vector<double> &getSsf() const;
// const Vector2D &getLfc() const;
// const std::vector<double> &getWvg() const;
// const Vector2D &getIdr() const;
// std::vector<double> getSdr() const;
// double getUInt() const;
// double getError() const;
```

Remove these from private section:

```cpp
// DELETE these lines:
// int runGrid() override;
// double getCoupling(GridPoint p) const override;
// double getDegeneracy(GridPoint p) const override;
// double getFxcIntegrandValue(GridPoint p) const override;
```

Add the `grid()` override:

```cpp
protected:
  // Grid accessor (override from VSBase)
  VSManager& grid() override { return grid_; }
  const VSManager& grid() const override { return grid_; }

private:
  std::shared_ptr<const VSStlsInput> inPtr;
  VSStlsManager grid_;  // ← renamed from 'grid' to 'grid_' to avoid name clash
```

**File**: `src/vs/vsstls.cpp`

1. Update constructor to use `grid_`:

```cpp
VSStls::VSStls(const std::shared_ptr<const VSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid_(in) {  // ← renamed
  if (in->getDimension() == dimensionsUtil::Dimension::D2) {
    throwError("2D calculations are not implemented for this scheme.");
  }
  setRsGrid();
  setFxcIntegrand();
}
```

2. Remove these method implementations (now in VSBase):
   - `runGrid()`
   - `getCoupling(GridPoint)`
   - `getDegeneracy(GridPoint)`
   - `getFxcIntegrandValue(GridPoint)`
   - `getSsf()`
   - `getLfc()`
   - `getWvg()`
   - `getIdr()`
   - `getSdr()`
   - `getUInt()`
   - `getError()`

3. Keep only:
   - `in()` implementation
   - `inScheme()` implementation
   - `computeQRaw(GridPoint)` implementation (which needs updating to use `grid_`)

Update `computeQRaw()` to use `grid_`:

```cpp
double VSStls::computeQRaw(GridPoint p) const {
  return QAdder::classical(
             grid_.VSManager::getWvg(p), grid_.VSManager::getSsf(p), inPtr)
      .get();
}
```

---

### 6. `VSQstlsManager` — Add scheme-specific convenience methods

**File**: `include/vs/qvsstls.hpp`

Add to `VSQstlsManager` public section:

```cpp
public:
  // Scheme-specific getters (downcast to Qstls)
  const Vector2D& getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
```

**File**: `src/vs/qvsstls.cpp`

Implement in `VSQstlsManager`:

```cpp
const Vector2D& VSQstlsManager::getIdr() const {
  return dynamic_cast<const Qstls&>(getWorkerAt(GridPoints::CENTER)).getIdr();
}

std::vector<double> VSQstlsManager::getSdr() const {
  return dynamic_cast<const Qstls&>(getWorkerAt(GridPoints::CENTER)).getSdr();
}

double VSQstlsManager::getUInt() const {
  return dynamic_cast<const Qstls&>(getWorkerAt(GridPoints::CENTER)).getUInt();
}
```

---

### 7. `QVSStls` — Implement `grid()` and remove duplicated methods

**File**: `include/vs/qvsstls.hpp`

Remove these method declarations from `QVSStls` (now inherited from `VSBase`):

```cpp
// DELETE these lines:
// const std::vector<double> &getSsf() const;
// const Vector2D &getLfc() const;
// const std::vector<double> &getWvg() const;
// const Vector2D &getIdr() const;
// std::vector<double> getSdr() const;
// double getUInt() const;
// double getError() const;
```

Remove these from private section:

```cpp
// DELETE these lines:
// int runGrid() override;
// double getCoupling(GridPoint p) const override;
// double getDegeneracy(GridPoint p) const override;
// double getFxcIntegrandValue(GridPoint p) const override;
```

Add the `grid()` override:

```cpp
protected:
  // Grid accessor (override from VSBase)
  VSManager& grid() override { return grid_; }
  const VSManager& grid() const override { return grid_; }

private:
  std::shared_ptr<const QVSStlsInput> inPtr;
  VSQstlsManager grid_;  // ← renamed from 'grid' to 'grid_' to avoid name clash
  std::shared_ptr<Integrator2D> itg2D;
  std::vector<double> itgGrid;
```

**File**: `src/vs/qvsstls.cpp`

1. Update constructor to use `grid_`:

```cpp
QVSStls::QVSStls(const std::shared_ptr<const QVSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid_(in),  // ← renamed
      itg2D(make_shared<Integrator2D>(ItgType::DEFAULT, in->getIntError())) {
  if (in->getDegeneracy() == 0.0) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  const bool segregatedItg = in->getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = grid_.VSManager::getWvg(GridPoints::CENTER); }
  setRsGrid();
  setFxcIntegrand();
}
```

2. Remove these method implementations (now in VSBase):
   - `runGrid()`
   - `getCoupling(GridPoint)`
   - `getDegeneracy(GridPoint)`
   - `getFxcIntegrandValue(GridPoint)`
   - `getSsf()`
   - `getLfc()`
   - `getWvg()`
   - `getIdr()`
   - `getSdr()`
   - `getUInt()`
   - `getError()`

3. Keep only:
   - `in()` implementation
   - `inScheme()` implementation
   - `computeQRaw(GridPoint)` implementation (which needs updating to use `grid_`)

Update `computeQRaw()` to use `grid_`:

```cpp
double QVSStls::computeQRaw(GridPoint p) const {
  const auto &w = dynamic_cast<const VSQstlsWorker &>(grid_.getWorkerAt(p));
  return w.computeQAdder(itg2D, itgGrid);
}
```

---

## Alternative Implementation Note

If the direct upcast `return grid_;` doesn't work due to C++ type compatibility issues (since `VSStlsManager` inherits from both `VSManager` and `Stls`), use:

```cpp
// In VSStls
VSManager& VSStls::grid() {
  return *StlsUtil::dynamic_pointer_cast<VSStlsManager, VSManager>(&grid_);
}
const VSManager& VSStls::grid() const {
  return *StlsUtil::dynamic_pointer_cast<const VSStlsManager, const VSManager>(&grid_);
}

// In QVSStls
VSManager& QVSStls::grid() {
  return *StlsUtil::dynamic_pointer_cast<VSQstlsManager, VSManager>(&grid_);
}
const VSManager& QVSStls::grid() const {
  return *StlsUtil::dynamic_pointer_cast<const VSQstlsManager, const VSManager>(&grid_);
}
```

---

## Summary of Files Changed

| File | Changes |
|------|---------|
| `include/vs/qvsstls.hpp` | Remove debug print from `VSQstlsWorker::init()`, add convenience methods to `VSQstlsManager`, remove duplicated methods from `QVSStls`, add `grid()` override, rename member `grid` → `grid_` |
| `include/vs/vsstls.hpp` | Add convenience methods to `VSStlsManager`, remove duplicated methods from `VSStls`, add `grid()` override, rename member `grid` → `grid_` |
| `include/vs/vsbase.hpp` | Add virtual `grid()` methods, remove virtual declarations for now-unified methods, add public getters |
| `include/vs/vsmanager.hpp` | Add convenience methods `getSsf()`, `getLfc()`, `getWvg()` |
| `src/vs/vsbase.cpp` | Implement all delegating methods using `grid()` |
| `src/vs/vsmanager.cpp` | Implement convenience methods |
| `src/vs/vsstls.cpp` | Remove duplicated method implementations, add convenience methods to `VSStlsManager`, update references to use `grid_` |
| `src/vs/qvsstls.cpp` | Remove duplicated method implementations, add convenience methods to `VSQstlsManager`, update references to use `grid_` |

---

## Benefits

1. **Eliminates code duplication** — 11 methods moved from both child classes to VSBase (22 → 11 implementations)
2. **Consistent with existing patterns** — mirrors the `in()` virtual method pattern
3. **No changes to VSWorker interface** — type-specific logic stays in manager classes
4. **Cleaner separation of concerns** — VSBase handles generic delegation, managers handle type-specific downcasting
5. **Easier to maintain** — changes to delegation logic only need to happen in one place

---

## Methods Eliminated from Child Classes

**From both `VSStls` and `QVSStls` (moved to `VSBase`):**
- `runGrid()`
- `getCoupling(GridPoint)`
- `getDegeneracy(GridPoint)`
- `getFxcIntegrandValue(GridPoint)`
- `getSsf()`
- `getLfc()`
- `getWvg()`
- `getIdr()`
- `getSdr()`
- `getUInt()`
- `getError()`

**Total reduction:** 22 method implementations → 11 in VSBase + 6 in managers (3 per manager class)
