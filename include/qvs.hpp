#ifndef QVS_HPP
#define QVS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstls.hpp"
#include "vector2D.hpp"
#include "vsbase.hpp"
#include <cmath>
#include <memory>

// // -----------------------------------------------------------------
// // Class to handle simultaneous state point calculations
// // -----------------------------------------------------------------

// class QStlsCSR : public CSR<Vector2D, Qstls, QVSStlsInput> {

// public:

//   // Constructor
//   explicit QStlsCSR(const QVSStlsInput &in_)
//       : CSR(in_, Qstls(in_.toQstlsInput(), false, false)) {}
//   // Compute auxiliary density response
//   void computeAdrStls();
//   void computeAdr();
//   // Update the static structure factor
//   void updateSsf() { ssf = ssfOld; };
//   // Initialize the scheme
//   void init();
//   // Compute Q
//   double getQAdder() const;

// private:

//   // Helper methods to compute the derivatives
//   double getDerivative(const std::shared_ptr<Vector2D> &f,
//                        const int &l,
//                        const size_t &idx,
//                        const Derivative &type);
// };

// // -----------------------------------------------------------------
// // Class to handle the Q-adder in the free parameter expression
// // -----------------------------------------------------------------

// class QAdder {

// public:

//   // Constructor
//   QAdder(const double &Theta_,
//          const double &mu_,
//          const double &limitMin,
//          const double &limitMax,
//          const std::vector<double> &itgGrid_,
//          Integrator1D &itg1_,
//          Integrator2D &itg2_,
//          const Interpolator1D &interp_)
//       : Theta(Theta_),
//         mu(mu_),
//         limits(limitMin, limitMax),
//         itgGrid(itgGrid_),
//         itg1(itg1_),
//         itg2(itg2_),
//         interp(interp_) {}
//   // Get Q-adder
//   double get() const;

// private:

//   const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
//   // Degeneracy parameter
//   const double Theta;
//   // Chemical potential
//   const double mu;
//   // Integration limits
//   const std::pair<double, double> limits;
//   // Grid for 2D integration
//   const std::vector<double> &itgGrid;
//   // Integrator objects
//   Integrator1D &itg1;
//   Integrator2D &itg2;
//   // Interpolator 1D class instance
//   const Interpolator1D &interp;

//   // SSF interpolation
//   double ssf(const double &y) const;
//   // Integrands
//   double integrandDenominator(const double q) const;
//   double integrandNumerator1(const double q) const;
//   double integrandNumerator2(const double w) const;
//   // Get Integral denominator
//   void getIntDenominator(double &res) const;
// };

// // -----------------------------------------------------------------
// // Class to handle the state point derivatives
// // -----------------------------------------------------------------

// class QStructProp : public StructPropBase<QStlsCSR, QVSStlsInput> {

// public:

//   // Constructor
//   explicit QStructProp(const QVSStlsInput &in)
//       : StructPropBase(in) {}
//   // Get Q term
//   std::vector<double> getQ() const;

// private:

//   using StructPropBase = StructPropBase<QStlsCSR, QVSStlsInput>;
//   // Setup input for the CSR objects
//   std::vector<QVSStlsInput> setupCSRInput(const QVSStlsInput &in);
//   // Setup dependencies in the CSR objects
//   void setupCSRDependencies();
//   // Perform iterations to compute structural properties
//   void doIterations();
// };

// // -----------------------------------------------------------------
// // Class to handle the thermodynamic properties and derivatives
// // -----------------------------------------------------------------

// class QThermoProp : public ThermoPropBase<QStructProp, QVSStlsInput> {

// public:

//   // Typdef
//   using ThermoPropBase = ThermoPropBase<QStructProp, QVSStlsInput>;
//   // Constructors
//   explicit QThermoProp(const QVSStlsInput &in)
//       : ThermoPropBase(in) {}
//   // Get internal energy and internal energy derivatives
//   std::vector<double> getQData() const;
//   // Get structural properties
//   Vector2D getAdr();
// };

// // -----------------------------------------------------------------
// // Solver for the qVS-STLS scheme
// // -----------------------------------------------------------------

// class QVSStls : public VSBase<QThermoProp, Rpa, QVSStlsInput> {

// public:

//   // Typedef
//   using VSBase = VSBase<QThermoProp, Rpa, QVSStlsInput>;
//   // Constructor from initial data
//   explicit QVSStls(const QVSStlsInput &in_)
//       : VSBase(in_) {}
//   // Constructor for recursive calculations
//   QVSStls(const QVSStlsInput &in_, const QThermoProp &thermoProp_)
//       : VSBase(in_, thermoProp_) {}
//   // Getters
//   Vector2D getAdr() const { return adr; }

// private:

//   // Compute free parameter
//   double computeAlpha();
//   // Iterations to solve the qvsstls-scheme
//   void updateSolution();
//   // Setup free energy integrand
//   void initFreeEnergyIntegrand();
//   // Auxiliary density response
//   Vector2D adr;
// };

// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class QVSStls: public VSBase {

public:

  // Constructor
  QVSStls(const QVSStlsInput& in_): VSBase(in_) {}
  // Getters
  Vector2D getAdr() const { return adr; }

private:

  // Compute free parameter
  double computeAlpha() { return 0.0; }
  // Iterations to solve the qvsstls-scheme
  void updateSolution() {};
  // Setup free energy integrand
  void initScheme() {};
  void initFreeEnergyIntegrand() {};
  // Auxiliary density response
  Vector2D adr;
};

#endif
