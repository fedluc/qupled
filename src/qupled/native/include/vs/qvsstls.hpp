#ifndef VS_QVSSTLS_HPP
#define VS_QVSSTLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstls.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"
#include <memory>

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

class VSQstlsWorker : public VSWorker, public Qstls {
public:

  using InputType = QVSStlsInput;

  VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in, GridPoint p);

  const Vector2D &getLfc() const override { return Qstls::getLfc(); }
  const std::vector<double> &getWvg() const override { return Qstls::getWvg(); }
  const std::vector<double> &getSsf() const override { return Qstls::getSsf(); }

  void init() override {
    std::cerr << "Initializing VSQstlsWorker..." << std::endl;
    Qstls::init();
  }
  void initialGuess() override { Qstls::initialGuess(); }
  void computeSsf() override { Qstls::computeSsf(); }
  void computeLfc() override { Qstls::computeLfc(); }
  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  double computeError() const override { return Qstls::computeError(); }
  void updateSolution() override { Qstls::updateSolution(); }

  double computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                       const std::vector<double> &itgGrid) const;
};

// -----------------------------------------------------------------
// VSQstlsManager
// -----------------------------------------------------------------

class VSQstlsManager : public VSManager, public Qstls {
public:

  explicit VSQstlsManager(const std::shared_ptr<const QVSStlsInput> &in);

  using VSManager::computeError;
  using VSManager::computeLfc;
  using VSManager::computeSsf;
  using VSManager::init;
  using VSManager::initialGuess;
  using VSManager::updateSolution;
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
  VSQstlsManager grid;
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
