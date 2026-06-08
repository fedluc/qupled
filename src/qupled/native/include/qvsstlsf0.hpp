#ifndef QVSSTLSF0_HPP
#define QVSSTLSF0_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstlsf0.hpp"
#include "vector2D.hpp"
#include "vsbase.hpp"
#include <memory>

class QThermoPropF0;
class QstlsF0CSR;
class QAdderF0;

class QVSStlsF0 : public VSBase, public QstlsF0 {

public:

  explicit QVSStlsF0(const std::shared_ptr<const QVSStlsF0Input> &in_);
  using VSBase::compute;

private:

  std::shared_ptr<QThermoPropF0> thermoProp;

  const VSInput &in() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
  }

  void init() override;
  double computeAlpha() override;
  void updateSolution() override;

  void print(const std::string &msg) { VSBase::print(msg); }
  void println(const std::string &msg) { VSBase::println(msg); }
};

class QThermoPropF0 : public ThermoPropBase {

public:

  explicit QThermoPropF0(const std::shared_ptr<const QVSStlsF0Input> &in_);
  std::vector<double> getQData() const;

private:

  std::shared_ptr<QstlsF0CSR> structProp;
};

class QstlsF0CSR : public CSR, public QstlsF0 {

public:

  explicit QstlsF0CSR(const std::shared_ptr<const QVSStlsF0Input> &in_)
      : QstlsF0CSR(in_, true) {}
  QstlsF0CSR(const std::shared_ptr<const QVSStlsF0Input> &in_, const bool isMaster_);

  int compute() override { return QstlsF0::compute(); }
  double getQAdder(const size_t &idx) const;

private:

  const std::shared_ptr<Integrator2D> itg2D;
  std::vector<double> itgGrid;

  const VSInput &inVS() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
  }

  const Input &inRpa() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, Input>(inPtr);
  }

  void setupWorkers(const QVSStlsF0Input &in) {
    CSR::setupWorkers<QstlsF0CSR, QVSStlsF0Input>(in);
  }

  void init() override { CSR::init(); }
  void computeLfc() override { CSR::computeLfc(); }
  void computeSsf() override { CSR::computeSsf(); }
  void initialGuess() override { CSR::initialGuess(); }
  void updateSolution() override { CSR::updateSolution(); }
  double computeError() const override { return CSR::computeError(); }

  void initWorker() override;
  void computeLfcWorker() override { QstlsF0::computeLfc(); }
  void computeSsfWorker() override { QstlsF0::computeSsf(); }
  void initialGuessWorker() override { QstlsF0::initialGuess(); }
  void updateSolutionWorker() override { QstlsF0::updateSolution(); }
  double computeErrorWorker() const override { return QstlsF0::computeError(); }
  Vector2D &getLfc() override { return lfc; }

  const std::vector<double> &getSsf() const override { return QstlsF0::getSsf(); }
  const std::vector<double> &getWvg() const override { return QstlsF0::getWvg(); }
  const Vector2D &getLfc() const override { return QstlsF0::getLfc(); }
};

class QAdderF0 {

public:

  QAdderF0(const double &rs_,
           const double &Theta_,
           const double &mu_,
           const double &etaIn,
           const double &ySecIn,
           const double &aCutoffIn,
           const double &limitMin,
           const double &limitMax,
           const std::vector<double> &itgGrid_,
           std::shared_ptr<Integrator1D> itg1_,
           std::shared_ptr<Integrator2D> itg2_,
           std::shared_ptr<Interpolator1D> interp_)
      : f0(rs_, Theta_, mu_, etaIn, ySecIn, aCutoffIn),
        limits(limitMin, limitMax),
        itgGrid(itgGrid_),
        itg1(itg1_),
        itg2(itg2_),
        interp(interp_) {}

  double get() const;

private:

  const QstlsF0Util::Distribution f0;
  const std::pair<double, double> limits;
  const std::vector<double> &itgGrid;
  const std::shared_ptr<Integrator1D> itg1;
  const std::shared_ptr<Integrator2D> itg2;
  const std::shared_ptr<Interpolator1D> interp;

  double ssf(const double &y) const;
  double integrandDenominator(const double q) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  void getIntDenominator(double &res) const;
};

#endif
