#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "input.hpp"
#include "state_point_grid.hpp"
#include "stls.hpp"
#include "vsbase.hpp"

// -----------------------------------------------------------------
// VSStlsWorker: Stls extended with VS-specific methods for StatePointGrid
// -----------------------------------------------------------------

class VSStlsWorker : public Stls {

public:

  using InputType = VSStlsInput;

  VSStlsWorker(const std::shared_ptr<const VSStlsInput> &in,
               bool verbose,
               GridPoint /* unused */)
      : Stls(in, verbose) {}

  // VS-specific: apply alpha-correction to lfc (called by StatePointGrid step
  // 3)
  void applyLfcDiff(const Vector2D &v) { lfc.diff(v); }

  // Expose protected Stls iteration steps for StatePointGrid orchestration
  void doInit() { Stls::init(); }
  void doInitialGuess() { Stls::initialGuess(); }
  void doComputeLfc() { Stls::computeLfc(); }
  void doComputeSsf() { Stls::computeSsf(); }
  void doUpdateSolution() { Stls::updateSolution(); }
  double doComputeError() const { return Stls::computeError(); }
};

// -----------------------------------------------------------------
// VSStls: VS scheme based on Stls structural properties
// -----------------------------------------------------------------

class VSStls : public VSBase {

public:

  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in);
  using VSBase::compute;

  // Public interface required by Python bindings
  const std::vector<double> &getSsf() const { return ssf; }
  const Vector2D &getLfc() const { return lfc; }
  const std::vector<double> &getWvg() const;
  const Vector2D &getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  double getError() const;

private:

  std::shared_ptr<const VSStlsInput> inPtr;
  StatePointGrid<VSStlsWorker> grid;
  std::vector<double> ssf;
  Vector2D lfc;

  const VSInput &in() const override;
  const Input &inScheme() const override;
  void init() override;
  void updateSolution() override;

  int runGrid() override;
  double getCoupling(GridPoint p) const override;
  double getDegeneracy(GridPoint p) const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  // Returns {uint, uintr, uintt} — classical limit of the Q-term
  std::vector<double> computeQData() override;
};

#endif
