#ifndef ESA_HPP
#define ESA_HPP

#include "auto_diff.hpp"
#include "rpa.hpp"

// -----------------------------------------------------------------
// Solver for the ESA scheme
// -----------------------------------------------------------------

class ESA : public Rpa {

public:

  // ESA constructor
  explicit ESA(const RpaInput &in_)
      : Rpa(in_) {}
  // Compute the scheme
  int compute();

private:

  // Static local field correction
  void computeSlfc();
  // Parametrization of the slfc obtained from neural networks
  double slfcNN(const double &x) const;
  // slfc from the compressibility sum rule
  double slfcCSR(const double &x) const;
  // Compute static local field correction coefficients
  struct SlfcCoefficients {
    // Coefficients for the long wavelength limit
    double lwl;
    // Coefficients for the activation function
    double afEta;
    double afxm;
    // Coefficients for the neural network parametrization
    double nna;
    double nnb;
    double nnc;
    double nnd;
    // Coefficients for the compressibility sum-rule
    double csr;
  };
  SlfcCoefficients slfcCoeff;
  void computeSlfcCoefficients();
  void computeSlfcNNCoefficients();
  void computeSlfcCSRCoefficients();
  // On top value of the radial distribution function
  double onTop() const;
  // Activation function for the asymptotic limit of slfc
  double activationFunction(const double &x) const;
  // Parametrization of the free energy
  AutoDiff2 freeEnergy(const double &rs, const double &theta) const;
  AutoDiff2 freeEnergy(const AutoDiff2 &rs, const AutoDiff2 &theta) const;
};

#endif
