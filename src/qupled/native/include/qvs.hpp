#ifndef QVS_HPP
#define QVS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstls.hpp"
#include "vector2D.hpp"
#include "vsbase.hpp"
#include <cmath>
#include <memory>

class QThermoProp;
class QStructProp;
class QstlsCSR;
class QAdder;

// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class QVSStls : public VSBase, public Qstls {

public:

  // Constructor from initial data
  explicit QVSStls(const QVSStlsInput &in_);
  // Constructor for recursive calculations
  QVSStls(const QVSStlsInput &in_, const QThermoProp &thermoProp_);
  // Solve the scheme
  using VSBase::compute;
  // Getters
  const QVSStlsInput &getInput() const { return in; }

private:

  // Input
  QVSStlsInput in;
  // Thermodyanmic properties
  std::shared_ptr<QThermoProp> thermoProp;
  // Initialize
  void initScheme();
  void initFreeEnergyIntegrand();
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the qvsstls-scheme
  void updateSolution();
  // Print info
  void print(const std::string &msg) { VSBase::print(msg); }
  void println(const std::string &msg) { VSBase::println(msg); }
};

// -----------------------------------------------------------------
// Class to handle the thermodynamic properties
// -----------------------------------------------------------------

class QThermoProp : public ThermoPropBase {

public:

  // Constructors
  explicit QThermoProp(const QVSStlsInput &in_);
  // Get internal energy and internal energy derivatives
  std::vector<double> getQData() const;
  // Get structural properties
  const Vector2D &getAdr();

private:

  std::shared_ptr<QStructProp> structProp;
};

// -----------------------------------------------------------------
// Class to handle the structural properties
// -----------------------------------------------------------------

class QStructProp : public Logger, public StructPropBase {

public:

  // Constructor
  explicit QStructProp(const QVSStlsInput &in_);
  // Get structural properties for output
  const QstlsCSR &getCsr(const Idx &idx) const;
  // Get Q term
  std::vector<double> getQ() const;

private:

  // Inputs
  const QVSStlsInput in;
  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<std::shared_ptr<QstlsCSR>> csr;
  // Setup dependencies in the CSR objects
  std::vector<QVSStlsInput> setupCSRInput();
  void setupCSR();
  // Perform iterations to compute structural properties
  void doIterations();
};

// -----------------------------------------------------------------
// Class to solve one state point
// -----------------------------------------------------------------

class QstlsCSR : public CSR, public Qstls {

public:

  // Constructor
  explicit QstlsCSR(const QVSStlsInput &in_);
  // Compute auxiliary density response
  void computeAdrQStls();
  void computeAdr();
  // Publicly esposed private stls methods
  void init();
  void initialGuess() { Qstls::initialGuess(); }
  void computeSsf() { Qstls::computeSsf(); }
  double computeError() { return Qstls::computeError(); }
  void updateSolution() { Qstls::updateSolution(); }
  // Update the static structure factor
  void updateSsf() { ssf = ssfOld; };
  // Compute Q
  double getQAdder() const;
  // Getters
  const std::vector<double> &getSsf() const { return Qstls::getSsf(); }
  const std::vector<double> &getSlfc() const { return Qstls::getSlfc(); }
  const std::vector<double> &getWvg() const { return Qstls::getWvg(); }

private:

  // Input
  QVSStlsInput in;
  // Compute derivative contribution to the auxiliary density response
  Vector2D getDerivativeContribution() const;
};

// -----------------------------------------------------------------
// Class to handle the Q-adder in the free parameter expression
// -----------------------------------------------------------------

class QAdder {

public:

  // Constructor
  QAdder(const double &Theta_,
         const double &mu_,
         const double &limitMin,
         const double &limitMax,
         const std::vector<double> &itgGrid_,
         Integrator1D &itg1_,
         Integrator2D &itg2_,
         const Interpolator1D &interp_)
      : Theta(Theta_),
        mu(mu_),
        limits(limitMin, limitMax),
        itgGrid(itgGrid_),
        itg1(itg1_),
        itg2(itg2_),
        interp(interp_) {}
  // Get Q-adder
  double get() const;

private:

  const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Integration limits
  const std::pair<double, double> limits;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  // Integrator objects
  Integrator1D &itg1;
  Integrator2D &itg2;
  // Interpolator 1D class instance
  const Interpolator1D &interp;

  // SSF interpolation
  double ssf(const double &y) const;
  // Integrands
  double integrandDenominator(const double q) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  // Get Integral denominator
  void getIntDenominator(double &res) const;
};

#endif
