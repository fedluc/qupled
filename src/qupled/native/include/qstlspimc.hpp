#ifndef QSTLSPIMC_HPP
#define QSTLSPIMC_HPP

#include "qstls.hpp"
#include <memory>

// -----------------------------------------------------------------
// Solver for the qSTLS-PIMC scheme
// -----------------------------------------------------------------

class QstlsPimc : public Qstls {

public:

  // Constructor
  explicit QstlsPimc(const std::shared_ptr<const QstlsPimcInput> &in_);

  // Getters
  const std::vector<double> &getF0Grid() const { return f0Grid; }
  const std::vector<double> &getF0Values() const { return f0Values; }

protected:

  // Initialize basic properties
  void init() override;
  // Compute static structure factor at finite temperature without HF split
  void computeSsfFinite() override;

private:

  // Input parameters
  const QstlsPimcInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, QstlsPimcInput>(inPtr);
  }
  // Stored interacting momentum distribution for plotting/output
  std::vector<double> f0Grid;
  std::vector<double> f0Values;
  // Compute fixed auxiliary density response
  void computeAdrFixedPimc();
};

namespace QstlsPimcUtil {

  // Interacting momentum distribution f0(y) and derivative from PIMC data
  class Distribution {

  public:

    Distribution(const double &rs,
                 const double &theta,
                 const double &etaIn,
                 const double &ySecIn,
                 const double &aCutoffIn);

    // Evaluate f0(y) and df0(y)/dy
    double value(const double &y) const;
    double derivative(const double &y) const;

  private:

    // State point
    const double rs;
    const double theta;
    // Tabulated PIMC data
    std::vector<double> yGrid;
    std::vector<double> fGrid;
    std::vector<double> dfGrid;
    // Domain of tabulated data
    double yMin;
    double yMax;
    double yJoin;
    // Switching function parameters
    double eta;
    double ySec;
    double aCutoff;
    // Asymptotic tail prefactor in coeff / y^8
    double tailCoeff;
    // Evaluate components
    double pimcValue(const double &y) const;
    double pimcDerivative(const double &y) const;
    double asymptoticValue(const double &y) const;
    double asymptoticDerivative(const double &y) const;
    double switchFunction(const double &y) const;
    double switchDerivative(const double &y) const;
    double determineYSec() const;
    // Moment normalization
    double value(const double &y, const double &etaLocal, const double &ySecLocal) const;
    double moment(const double &etaLocal, const double &ySecLocal) const;
    double solveYSec(const double &etaLocal) const;
    double solveEta(const double &ySecLocal) const;
  };

  class AdrFixedPimc : public QstlsUtil::AdrFixedBase {

  public:

    // Constructor for finite temperature calculations
    AdrFixedPimc(const double &Theta_,
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

    // Get integration result
    void get(const std::vector<double> &wvg, Vector3D &res) const;

  private:

    // Integrands
    double integrand1(const double &q, const double &l) const;
    double integrand2(const double &t, const double &y, const double &l) const;
    // Interacting momentum distribution
    const Distribution &f0;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
  };

} // namespace QstlsPimcUtil

#endif
