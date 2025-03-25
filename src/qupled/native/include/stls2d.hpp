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
  void computeSsf2D();
  // Iterations to solve the stls scheme
  void doIterations();
};

// -----------------------------------------------------------------
// Class to handle the 2D HF SSF 
// -----------------------------------------------------------------

class SSFHF2D {

  public:
  
    // Constructor
    SSFHF2D(const double &Theta_,
           const double &x_,
           const double &mu_,
           const double &limitMin,
           const double &limitMax,
           const std::vector<double> &itgGrid_,
           Integrator2D &itg2_)
        : Theta(Theta_),
          x(x_),
          mu(mu_),
          limits(limitMin, limitMax),
          itgGrid(itgGrid_),
          itg2(itg2_) {}
    // Get Q-adder
    double get() const;
  
  private:
  
    // Degeneracy parameter
    const double Theta;
    const double x;
    // Chemical potential
    const double mu;
    // Integration limits
    const std::pair<double, double> limits;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
    // Integrator objects
    Integrator2D &itg2;

    // Integrands
    double integrandOut(const double q) const;
    double integrandIn(const double w) const;
  };

#endif