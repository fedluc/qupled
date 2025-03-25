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

class Stls2D : public Rpa {

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

#endif