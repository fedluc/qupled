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
  Rpa(const std::shared_ptr<const Input> &in_, const bool verbose_);
  explicit Rpa(const std::shared_ptr<const Input> &in_)
      : Rpa(in_, true) {}

protected:

  // Hartree-Fock Static structure factor
  std::vector<double> ssfHF;
  // Initialize basic properties
  void init() override;
  // Compute static structure factor
  void computeSsfFinite() override;
  void computeSsfGround() override;

private:

  // Compute Hartree-Fock static structure factor
  void computeSsfHF();
  // Compute local field correction
  void computeLfc() override;

};

namespace RpaUtil {

  // -----------------------------------------------------------------
  // Classes for the static structure factor
  // -----------------------------------------------------------------

  class SsfBase {

  protected:

    // Constructor
    SsfBase(const double &x_,
            const double &Theta_,
            const double &rs_,
            const double &ssfHF_,
            std::span<const double> lfc_)
        : x(x_),
          Theta(Theta_),
          rs(rs_),
          ssfHF(ssfHF_),
          lfc(lfc_) {}
    // Wave-vector
    const double x;
    // Degeneracy parameter
    const double Theta;
    // Coupling parameter
    const double rs;
    // Hartree-Fock contribution
    const double ssfHF;
    // Local field correction
    std::span<const double> lfc;
    // Normalized interaction potential
    const double ip = 4.0 * numUtil::lambda * rs / (M_PI * x * x);
    // Normalized 2D interaction potential
    const double ip2D = sqrt(2.0) * rs / x;
  };

  class Ssf : public SsfBase, dimensionsUtil::DimensionsHandler {

  public:

    // Constructor
    Ssf(const double &x_,
        const double &Theta_,
        const double &rs_,
        const double &ssfHF_,
        std::span<const double> lfc_,
        std::span<const double> idr_,
        const std::shared_ptr<const Input> in_)
        : SsfBase(x_, Theta_, rs_, ssfHF_, lfc_),
          idr(idr_), in(in_), res(x_) {}
    // Get static structore factor
    double get();

  protected:

    // Ideal density response
    const std::span<const double> idr;
    const std::shared_ptr<const Input> in;

  private:
    
    double res;
    void compute2D() override;
    void compute3D() override;

  };

  class SsfGround : public SsfBase {

  public:

    // Constructor for zero temperature calculations
    SsfGround(const double &x_,
              const double &rs_,
              const double &ssfHF_,
              std::span<const double> lfc_,
              const double &OmegaMax_,
              std::shared_ptr<Integrator1D> itg_)
        : SsfBase(x_, 0, rs_, ssfHF_, lfc_),
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
