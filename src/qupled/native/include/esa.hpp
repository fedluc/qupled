#ifndef ESA_HPP
#define ESA_HPP

#include "rpa.hpp"
#include <cmath>

// Forward declarations
class Dual22;

// -----------------------------------------------------------------
// Solver for the ESA scheme
// -----------------------------------------------------------------

class ESA : public Rpa {

public:

  // ESA constructor
  explicit ESA(const Input &in_)
      : Rpa(in_) {}
  // Compute the scheme
  int compute();

private:

  // Static local field correction
  void computeSlfc();
};

namespace ESAUtil {

  class Slfc {

  public:

    // Constructor
    explicit Slfc(const double &rs_, const double &theta_)
        : rs(rs_),
          theta(theta_) {}
    // Get the static local field correction for a given wave-vector x
    double get(const double &x);

  public:

    // Constant for unit conversion
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
    // Coupling parameter
    const double rs;
    // Degeneracy parameter
    const double theta;
    // Compute static local field correction coefficients
    struct Coefficients {
      // Flag marking whether the coefficients are valid or should be recomputed
      bool valid = false;
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
    // Coefficients
    Coefficients coeff;
    // Parametrization of the slfc obtained from neural networks
    double nn(const double &x) const;
    // slfc from the compressibility sum rule
    double csr(const double &x) const;
    // Compute coefficients
    void computeCoefficients();
    void computeNNCoefficients();
    void computeCSRCoefficients();
    // On top value of the radial distribution function
    double onTop() const;
    // Activation function for the asymptotic limit of slfc
    double activationFunction(const double &x) const;
    // Parametrization of the free energy
    Dual22 freeEnergy() const;
    Dual22 freeEnergy(const Dual22 &rs, const Dual22 &theta) const;
  };

} // namespace ESAUtil
#endif
