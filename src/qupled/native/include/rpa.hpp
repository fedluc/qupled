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
            const double &ssfHF_,
            std::span<const double> lfc_,
            const std::shared_ptr<const Input> in_)
        : x(x_),
          ssfHF(ssfHF_),
          lfc(lfc_),
          in(in_) {}
    // Wave-vector
    const double x;
    // Hartree-Fock contribution
    const double ssfHF;
    // Local field correction
    std::span<const double> lfc;
    // Input struct
    const std::shared_ptr<const Input> in;
    // Normalized interaction potential
    double ip() const;
  };

  class Ssf : public SsfBase, dimensionsUtil::DimensionsHandler {

  public:

    // Constructor
    Ssf(const double &x_,
        const double &ssfHF_,
        std::span<const double> lfc_,
        const std::shared_ptr<const Input> in_,
        std::span<const double> idr_)
        : SsfBase(x_, ssfHF_, lfc_, in_),
          idr(idr_),
          res(numUtil::NaN) {}
    // Get static structore factor
    double get();

  protected:

    // Ideal density response
    const std::span<const double> idr;

  private:

    // Result of the integration
    double res;
    // Compute methods
    void compute2D() override;
    void compute3D() override;
  };

  class SsfGround : public SsfBase {

  public:

    // Constructor for zero temperature calculations
    SsfGround(const double &x_,
              const double &ssfHF_,
              std::span<const double> lfc_,
              std::shared_ptr<Integrator1D> itg_,
              const std::shared_ptr<const Input> in_)
        : SsfBase(x_, ssfHF_, lfc_, in_),
          itg(itg_) {}
    // Get result of integration
    double get();

  protected:

    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Integrand for zero temperature calculations
    double integrand(const double &Omega) const;
  };

} // namespace RpaUtil

#endif
