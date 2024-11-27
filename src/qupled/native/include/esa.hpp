#ifndef ESA_HPP
#define ESA_HPP

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

protected:

  // Static local field correction
  void computeSlfc();
  // On top value of the radial distribution function
  double onTop() const;
  // Activation function for the asymptotic limit of slfc
  double activationFunction(const double &x) const;
  // Parametrization of the slfc obtained from neural networks
  double slfcNN(const double &x) const;
  // slfc from the compressibility sum rule
  double slfcCSR(const double &x) const;
  // Parametrization of the free energy
  AutoDiff2D freeEnergy(const double &rs, const double &theta) const;
  AutoDiff2D freeEnergy(const AutoDiff2D &rs, const AutoDiff2D &theta) const;
};

#endif
