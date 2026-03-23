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
};

// -----------------------------------------------------------------
// VSStls: VS scheme based on Stls structural properties
// -----------------------------------------------------------------

class VSStls : public VSBase {
public:

  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in);
  using VSBase::compute;

  const std::vector<double> &getSsf() const;
  const Vector2D &getLfc() const;
  const std::vector<double> &getWvg() const;
  const Vector2D &getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  double getError() const;

private:

  std::shared_ptr<const VSStlsInput> inPtr;
  VSStlsManager grid;

  const VSInput &in() const override;
  const Input &inScheme() const override;

  int runGrid() override;
  double getCoupling(GridPoint p) const override;
  double getDegeneracy(GridPoint p) const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  double computeQRaw(GridPoint p) const override;
};

#endif
