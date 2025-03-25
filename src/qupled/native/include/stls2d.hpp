#ifndef STLS2D_HPP
#define STLS2D_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "rpa.hpp"
#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Solver for the STLS scheme
// -----------------------------------------------------------------

class Stls2D : public Stls {

public:

  // Constructors
  Stls2D(const StlsInput &in_, const bool verbose_, const bool writeFiles_);
  explicit Stls2D(const StlsInput &in_)
      : Stls2D(in_, true, true) {}
  // Compute stls scheme
  int compute();
  // Getters
  const StlsInput &getInput() const { return in; }
  double getError() const { return computeError(); }

protected:

  // Input parameters
  StlsInput in;
  // Flag to write the recovery files
  const bool writeFiles;
  // Static local field correction to use during the iterations
  std::vector<double> slfcNew;
  // Initialize basic properties
  void init();
  // Compute static local field correction
  void computeSlfc2D();
  void computeSlfcStls2D();
  void computeChemicalPotential2D();
  void computeIdr2D();
  void computeSsfHF2D();
  // Iterations to solve the stls scheme
  void doIterations();
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

#endif