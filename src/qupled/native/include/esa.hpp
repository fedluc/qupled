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

  // Funtion for the ESA static local field correction
  void computeSlfc();
  // Function for free energy derivatives
  double fxc(const double &theta, const double &rs) const;
  // Resolution for the free energy derivatives
  const double dx = 1e-3;
};

#endif
