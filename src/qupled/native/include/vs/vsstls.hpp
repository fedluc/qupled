#ifndef VS_VSSTLS_HPP
#define VS_VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmaster_base.hpp"

// -----------------------------------------------------------------
// VSStlsWorker: implements IVSWorker for STLS-based VS workers
// -----------------------------------------------------------------

class VSStlsWorker : public VSWorkerBase, public Stls {
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
// VSStlsMaster: manages 9 VSStlsWorkers, inherits Stls loop
// -----------------------------------------------------------------

class VSStlsMaster : public VSMasterBase, public Stls {
public:

  explicit VSStlsMaster(const std::shared_ptr<const VSStlsInput> &in);

protected:

  void init() override { VSMasterBase::init(); }
  void computeLfc() override { VSMasterBase::computeLfc(); }
  void computeSsf() override { VSMasterBase::computeSsf(); }
  double computeError() const override { return VSMasterBase::computeError(); }
  void updateSolution() override { VSMasterBase::updateSolution(); }
  void initialGuess() override { VSMasterBase::initialGuess(); }
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
  VSStlsMaster grid;

  const VSInput &in() const override;
  const Input &inScheme() const override;

  int runGrid() override;
  double getCoupling(GridPoint p) const override;
  double getDegeneracy(GridPoint p) const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  double computeQRaw(GridPoint p) const override;
};

#endif
