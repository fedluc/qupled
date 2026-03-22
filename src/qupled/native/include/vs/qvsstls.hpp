#ifndef VS_QVSSTLS_HPP
#define VS_QVSSTLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstls.hpp"
#include "vs/vs_master.hpp"
#include "vs/vsbase.hpp"
#include <memory>

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

class VSQstlsWorker : public VSWorkerBase, public Qstls {
public:

  using InputType = QVSStlsInput;

  VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in,
                bool verbose,
                GridPoint p);

  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  void computeBaseLfc() override { Qstls::computeLfc(); }

  const Vector2D &getLfc() const override { return lfc; }
  const std::vector<double> &getWvg() const override { return wvg; }
  const std::vector<double> &getSsf() const override { return ssf; }

  void workerInit() override { Qstls::init(); }
  void workerInitialGuess() override { Stls::initialGuess(); }
  void workerComputeSsf() override { Stls::computeSsf(); }
  double workerComputeError() const override { return Stls::computeError(); }
  void workerUpdateSolution() override { Stls::updateSolution(); }

  double computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                       const std::vector<double> &itgGrid) const;
};

// -----------------------------------------------------------------
// VSQstlsMaster
// -----------------------------------------------------------------

class VSQstlsMaster : public VSMasterBase, public Qstls {
public:

  explicit VSQstlsMaster(const std::shared_ptr<const QVSStlsInput> &in);

protected:

  void init() override { masterInit(); }
  void computeLfc() override { masterComputeLfc(); }
  void computeSsf() override { masterComputeSsf(); }
  double computeError() const override { return masterComputeError(); }
  void updateSolution() override { masterUpdateSolution(); }
  void initialGuess() override { masterInitialGuess(); }
};

// -----------------------------------------------------------------
// QVSStls
// -----------------------------------------------------------------

class QVSStls : public VSBase {
public:

  explicit QVSStls(const std::shared_ptr<const QVSStlsInput> &in);
  using VSBase::compute;

  const std::vector<double> &getSsf() const;
  const Vector2D &getLfc() const;
  const std::vector<double> &getWvg() const;
  const Vector2D &getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  double getError() const;

private:

  std::shared_ptr<const QVSStlsInput> inPtr;
  VSQstlsMaster grid;
  std::shared_ptr<Integrator2D> itg2D;
  std::vector<double> itgGrid;

  const VSInput &in() const override;
  const Input &inScheme() const override;

  int runGrid() override;
  double getCoupling(GridPoint p) const override;
  double getDegeneracy(GridPoint p) const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  double computeQRaw(GridPoint p) const override;
};

#endif
