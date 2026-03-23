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
  const Vector2D &getIdr() const override { return Qstls::getIdr(); }
  std::vector<double> getSdr() const override { return Qstls::getSdr(); }
  double getUInt() const override { return Qstls::getUInt(); }

  void init() override { Qstls::init(); }
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

  // Override from VSManager
  int compute() override { return Qstls::compute(); }
  double computeQRaw(GridPoint p) const override;

private:
  std::shared_ptr<Integrator2D> itg2D;
  std::vector<double> itgGrid;
};

// -----------------------------------------------------------------
// QVSStls
// -----------------------------------------------------------------

class QVSStls : public VSBase {
public:

  explicit QVSStls(const std::shared_ptr<const QVSStlsInput> &in);
  using VSBase::compute;

protected:

  VSManager &grid() override { return grid_; }
  const VSManager &grid() const override { return grid_; }

private:

  std::shared_ptr<const QVSStlsInput> inPtr;
  VSQstlsManager grid_;

  const VSInput &in() const override;
  const Input &inScheme() const override;
  double computeQRaw(GridPoint p) const override;
};

#endif
