#ifndef ESA_HPP
#define ESA_HPP

#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the ESA scheme
// -----------------------------------------------------------------

class ESA : public Stls {

private:

  // Function to perform the ESA scheme
  void doESA();
  // Funtion for the ESA static local field correction
  void computeSlfc();
  // Function for free energy derivatives
  double fxc(double theta, double rs) const;
  // Assign dx for the free energy derivatives
  const double dx = 1e-6; 

public:

  // ESA constructor
  ESA(const StlsInput& input) : Stls(input) { ; }
  // Compute ESA scheme
  int compute();
  
};

#endif