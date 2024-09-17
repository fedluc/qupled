#ifndef STLS_HPP
#define STLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "rpa.hpp"
#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Solver for the STLS scheme
// -----------------------------------------------------------------

class Stls : public Rpa {

public:

  // Constructors
  Stls(const StlsInput &in_, const bool verbose_, const bool writeFiles_);
  explicit Stls(const StlsInput &in_)
      : Stls(in_, true, true) {}
  // Compute stls scheme
  int compute();
  // Getters
  const StlsInput &getInput() const { return in; }
  double getError() const { return computeError(); }
  const std::vector<double> &getBf() const { return bf; }

protected:

  // Input parameters
  StlsInput in;
  // Flag to write the recovery files
  const bool writeFiles;
  // iet schemes
  bool useIet;
  // Static local field correction to use during the iterations
  std::vector<double> slfcNew;
  // Bridge function (for iet schemes)
  std::vector<double> bf;
  // Initialize basic properties
  void init();
  // Compute static local field correction
  void computeSlfc();
  void computeSlfcStls();
  void computeSlfcIet();
  // Compute bridge function
  void computeBf();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  bool initialGuessFromRecovery();
  bool initialGuessFromInput();
  double computeError() const;
  void updateSolution();
  // Write recovery files
  void writeRecovery();
  void readRecovery(std::vector<double> &wvgFile,
                    std::vector<double> &slfcFile) const;
};

// -----------------------------------------------------------------
// Classes for the static local field correction
// -----------------------------------------------------------------

class SlfcBase {

protected:

  // Wave-vector
  const double x;
  // Integration limits
  const double yMin;
  const double yMax;
  // Static structure factor interpolator
  const Interpolator1D &ssfi;
  // Compute static structure factor
  double ssf(const double &y) const;
  // Constructor
  SlfcBase(const double &x_,
           const double &yMin_,
           const double &yMax_,
           const Interpolator1D &ssfi_)
      : x(x_),
        yMin(yMin_),
        yMax(yMax_),
        ssfi(ssfi_) {}
};

class Slfc : public SlfcBase {

public:

  // Constructor
  Slfc(const double &x_,
       const double &yMin_,
       const double &yMax_,
       const Interpolator1D &ssfi_,
       Integrator1D &itg_)
      : SlfcBase(x_, yMin_, yMax_, ssfi_),
        itg(itg_) {}
  // Get result of integration
  double get() const;

private:

  // Integrator object
  Integrator1D &itg;
  // Integrand
  double integrand(const double &y) const;
};

class SlfcIet : public SlfcBase {

public:

  // Constructor
  SlfcIet(const double &x_,
          const double &yMin_,
          const double &yMax_,
          const Interpolator1D &ssfi_,
          const Interpolator1D &slfci_,
          const Interpolator1D &bfi_,
          const std::vector<double> &itgGrid_,
          Integrator2D &itg_)
      : SlfcBase(x_, yMin_, yMax_, ssfi_),
        itg(itg_),
        itgGrid(itgGrid_),
        slfci(slfci_),
        bfi(bfi_) {}
  // Get result of integration
  double get() const;

private:

  // Integrator object
  Integrator2D &itg;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  // Integrands
  double integrand1(const double &y) const;
  double integrand2(const double &w) const;
  // Static local field correction interpolator
  const Interpolator1D &slfci;
  // Bridge function interpolator
  const Interpolator1D &bfi;
  // Compute static local field correction
  double slfc(const double &x) const;
  // Compute bridge function
  double bf(const double &x_) const;
};

class BridgeFunction {

public:

  // Constructor
  BridgeFunction(const std::string &theory_,
                 const std::string &mapping_,
                 const double &rs_,
                 const double &Theta_,
                 const double &x_,
                 Integrator1D &itg_)
      : theory(theory_),
        mapping(mapping_),
        rs(rs_),
        Theta(Theta_),
        x(x_),
        itg(itg_) {}
  // Get result of the integration
  double get() const;

private:

  // Theory to be solved
  const std::string theory;
  // Iet mapping
  const std::string mapping;
  // Coupling parameter
  const double rs;
  // Degeneracy parameter
  const double Theta;
  // Wave vector
  const double x;
  // Integrator object
  Integrator1D &itg;
  // Constant for unit conversion
  const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
  // Hypernetted-chain bridge function
  double hnc() const;
  // Ichimaru bridge function
  double ioi() const;
  // Lucco Castello and Tolias bridge function
  double lct() const;
  double lctIntegrand(const double &r, const double &Gamma) const;
  // Coupling parameter to compute the bridge function
  double couplingParameter() const;
};

#endif
