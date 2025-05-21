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
  explicit ESA(const std::shared_ptr<const Input> &in_)
      : Rpa(in_) {}

private:

  // Local field correction
  void computeLfc() override;
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
    // Parametrization of the lfc obtained from neural networks
    double nn(const double &x) const;
    // lfc from the compressibility sum rule
    double csr(const double &x) const;
    // Compute coefficients
    void computeCoefficients();
    void computeNNCoefficients();
    void computeCSRCoefficients();
    // On top value of the radial distribution function
    double onTop() const;
    // Activation function for the asymptotic limit of lfc
    double activationFunction(const double &x) const;
    // Parametrization of the free energy
    Dual22 freeEnergy() const;
    Dual22 freeEnergy(const Dual22 &rs, const Dual22 &theta) const;
  };

} // namespace ESAUtil
#endif
