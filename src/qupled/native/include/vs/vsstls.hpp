#ifndef VS_VSSTLS_HPP
#define VS_VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"

// -----------------------------------------------------------------
// VSStlsWorker: implements VSWorker for STLS-based VS workers
// -----------------------------------------------------------------

class VSStlsWorker : public VSWorker, public Stls {
public:

  VSStlsWorker(const std::shared_ptr<const VSStlsInput> &in)
      : Stls(in, false) {}
  const Vector2D &getLfc() const override { return Stls::getLfc(); }
  const std::vector<double> &getWvg() const override { return Stls::getWvg(); }
  const std::vector<double> &getSsf() const override { return Stls::getSsf(); }
  const Vector2D &getIdr() const override { return Stls::getIdr(); }
  std::vector<double> getSdr() const override { return Stls::getSdr(); }
  double getUInt() const override { return Stls::getUInt(); }
  double getQAdder() const override { return getUInt(); }
  double getCoupling() const override { return inPtr->getCoupling(); }
  double getDegeneracy() const override { return inPtr->getDegeneracy(); }
  void init() override { Stls::init(); }
  void initialGuess() override { Stls::initialGuess(); }
  void computeSsf() override { Stls::computeSsf(); }
  void computeLfc() override { Stls::computeLfc(); }
  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  double computeError() const override { return Stls::computeError(); }
  void updateSolution() override { Stls::updateSolution(); }
};

// -----------------------------------------------------------------
// VSStlsManager: manages 9 VSStlsWorkers, inherits Stls loop
// -----------------------------------------------------------------

class VSStlsManager : public VSManager, public Stls {
public:

  explicit VSStlsManager(const std::shared_ptr<const VSStlsInput> &in);

  using VSManager::computeError;
  using VSManager::computeLfc;
  using VSManager::computeSsf;
  using VSManager::init;
  using VSManager::initialGuess;
  using VSManager::updateSolution;

  // Override from VSManager
  int compute() override { return Stls::compute(); }
  double computeQRaw(GridPoint p) const;

private:

  std::shared_ptr<const VSStlsInput> managerInPtr;
};

// -----------------------------------------------------------------
// VSStls: VS scheme based on Stls structural properties
// -----------------------------------------------------------------

class VSStls : public VSBase {
public:

  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in);
  using VSBase::compute;

private:

  std::shared_ptr<const VSStlsInput> inPtr;
  VSStlsManager grid_;

  VSManager &grid() override { return grid_; }
  const VSManager &grid() const override { return grid_; }
  const VSInput &in() const override;
  const Input &inScheme() const override;
};

#endif
