#ifndef ESA_HPP
#define ESA_HPP

#include "rpa.hpp"

// -----------------------------------------------------------------
// Solver for the ESA scheme
// -----------------------------------------------------------------

class ESA : public Rpa {
  
private:

  // Funtion for the ESA static local field correction
  void computeSlfc();
  // Function for free energy derivatives
  double fxc(double theta, double rs) const;
  // Resolution for the free energy derivatives
  const double dx = 1e-6; 

public:

  // ESA constructor
  ESA(const RpaInput& input);
  
};

#endif
