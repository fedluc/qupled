#ifndef QSTLSF0_HPP
#define QSTLSF0_HPP

#include "qstls.hpp"
#include <memory>

class QstlsF0 : public Qstls {

public:

  explicit QstlsF0(const std::shared_ptr<const QstlsPimcInput> &in_);

  const std::vector<double> &getF0Grid() const { return f0Grid; }
  const std::vector<double> &getF0Values() const { return f0Values; }

protected:

  void init() override;
  void computeSsfFinite() override;

  const QstlsPimcInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, QstlsPimcInput>(inPtr);
  }

private:

  std::vector<double> f0Grid;
  std::vector<double> f0Values;

  void computeAdrFixedF0();
};

namespace QstlsF0Util {

  class Distribution {

  public:

    Distribution(const double &rs,
                 const double &theta,
                 const double &mu,
                 const double &etaIn,
                 const double &ySecIn,
                 const double &aCutoffIn);

    double value(const double &y) const;
    double derivative(const double &y) const;

  private:

    const double rs;
    const double theta;
    const double mu;

    double eta;
    double ySec;
    double aCutoff;
    double tailCoeff;

    double fermiDiracValue(const double &y) const;
    double fermiDiracDerivative(const double &y) const;
    double asymptoticValue(const double &y) const;
    double asymptoticDerivative(const double &y) const;
    double switchFunction(const double &y) const;
    double switchDerivative(const double &y) const;
    double determineYSec() const;

    double value(const double &y, const double &etaLocal, const double &ySecLocal) const;
    double moment(const double &etaLocal, const double &ySecLocal) const;
    double solveYSec(const double &etaLocal) const;
  };

  class AdrFixedF0 : public QstlsUtil::AdrFixedBase {

  public:

    AdrFixedF0(const double &Theta_,
               const double &qMin_,
               const double &qMax_,
               const double &x_,
               const double &mu_,
               const Distribution &f0_,
               const std::vector<double> &itgGrid_,
               std::shared_ptr<Integrator2D> itg_)
        : QstlsUtil::AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
          f0(f0_),
          itg(itg_),
          itgGrid(itgGrid_) {}

    void get(const std::vector<double> &wvg, Vector3D &res) const;

  private:

    double integrand1(const double &q, const double &l) const;
    double integrand2(const double &t, const double &y, const double &l) const;

    const Distribution &f0;
    const std::shared_ptr<Integrator2D> itg;
    const std::vector<double> &itgGrid;
  };

} // namespace QstlsF0Util

#endif
