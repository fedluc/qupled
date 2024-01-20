#ifndef ESA_HPP
#define ESA_HPP

#include <limits>
#include <map>
#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the ESA scheme
// -----------------------------------------------------------------

class ESA : public Stls {

protected:

  // Function to perform the ESA scheme
  void doESA();

public:

  // ESA constructor (inherited from Stls)
  ESA(const StlsInput& input);
  // Compute ESA scheme
  int compute();
  // Funtion for the ESA SLFC
  void computeSlfc();
  // Derivatives of the QMC fxc
  double derivative_wrt_rs(double theta, double rs, double dx) const;
  double derivative_wrt_theta(double theta, double rs, double dx) const;
  double second_derivative_wrt_rs(double theta, double rs, double dx) const;
  double second_derivative_wrt_theta(double theta, double rs, double dx) const;
  double mixed_derivative(double theta, double rs, double dx) const;


};

#endif