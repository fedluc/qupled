#ifndef RPA_HPP
#define RPA_HPP

#include "hf.hpp"
#include "input.hpp"
#include "logger.hpp"
#include "numerics.hpp"
#include "vector2D.hpp"
#include <vector>

// -----------------------------------------------------------------
// Solver for the Random-Phase approximation scheme
// -----------------------------------------------------------------

class Rpa : public HF {

public:

  // Constructor
  Rpa(const Input &in_, const bool verbose_);
  explicit Rpa(const Input &in_)
      : Rpa(in_, true) {}
  // Compute the scheme
  int compute();

protected:

  // Hartree-Fock Static structure factor
  std::vector<double> ssfHF;
  // Initialize basic properties
  void init();
  // Compute Hartree-Fock static structure factor
  void computeSsfHF();
  // Compute static structure factor at finite temperature
  void computeSsf();
  void computeSsfFinite();
  void computeSsfGround();
  // Compute static local field correction
  void computeSlfc();
};

namespace RpaUtil {

  // -----------------------------------------------------------------
  // Classes for the static structure factor
  // -----------------------------------------------------------------

  class SsfBase {

  protected:

    // Wave-vector
    const double x;
    // Degeneracy parameter
    const double Theta;
    // Coupling parameter
    const double rs;
    // Hartree-Fock contribution
    const double ssfHF;
    // Static local field correction
    const double slfc;
    // Constant for unit conversion
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
    // Normalized interaction potential
    const double ip = 4.0 * lambda * rs / (M_PI * x * x);
    // Constructor
    SsfBase(const double &x_,
            const double &Theta_,
            const double &rs_,
            const double &ssfHF_,
            const double &slfc_)
        : x(x_),
          Theta(Theta_),
          rs(rs_),
          ssfHF(ssfHF_),
          slfc(slfc_) {}
  };

  class Ssf : public SsfBase {

  public:

    // Constructor
    Ssf(const double &x_,
        const double &Theta_,
        const double &rs_,
        const double &ssfHF_,
        const double &slfc_,
        const int nl_,
        const double *idr_)
        : SsfBase(x_, Theta_, rs_, ssfHF_, slfc_),
          nl(nl_),
          idr(idr_) {}
    // Get static structore factor
    double get() const;

  protected:

    // Number of Matsubara frequencies
    const int nl;
    // Ideal density response
    const double *idr;
  };

  class SsfGround : public SsfBase {

  public:

    // Constructor for zero temperature calculations
    SsfGround(const double &x_,
              const double &rs_,
              const double &ssfHF_,
              const double &slfc_,
              const double &OmegaMax_,
              std::shared_ptr<Integrator1D> itg_)
        : SsfBase(x_, 0, rs_, ssfHF_, slfc_),
          OmegaMax(OmegaMax_),
          itg(itg_) {}
    // Get result of integration
    double get();

  protected:

    // Integration limit
    const double OmegaMax;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Integrand for zero temperature calculations
    double integrand(const double &Omega) const;
  };

} // namespace RpaUtil

#endif
