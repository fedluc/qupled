#ifndef VS_VSSTLS_HPP
#define VS_VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vs/vs_master_base.hpp"
#include "vs/vsbase.hpp"

// -----------------------------------------------------------------
// VSStlsWorker: implements IVSWorker for STLS-based VS workers
// -----------------------------------------------------------------

class VSStlsWorker : public VSWorkerBase, public Stls {
public:

  using InputType = VSStlsInput;

  VSStlsWorker(const std::shared_ptr<const VSStlsInput> &in,
               bool verbose,
               GridPoint /* unused */)
      : Stls(in, verbose) {}

  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  void computeBaseLfc() override { Stls::computeLfc(); }

  const Vector2D &getLfc() const override { return lfc; }
  const std::vector<double> &getWvg() const override { return wvg; }
  const std::vector<double> &getSsf() const override { return ssf; }

  void init() override { HF::init(); }
  void initialGuess() override { Stls::initialGuess(); }
  void computeSsf() override { Stls::computeSsf(); }
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

  void init() override { masterInit(); }
  void computeLfc() override { masterComputeLfc(); }
  void computeSsf() override { masterComputeSsf(); }
  double computeError() const override { return masterComputeError(); }
  void updateSolution() override { masterUpdateSolution(); }
  void initialGuess() override { masterInitialGuess(); }
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
