#ifndef HF_HPP
#define HF_HPP

#include "dimensions_util.hpp"
#include "input.hpp"
#include "logger.hpp"
#include "numerics.hpp"
#include <vector>

// -----------------------------------------------------------------
// Solver for the Hartree-Fock scheme
// -----------------------------------------------------------------

class HF : public Logger {

public:

  // Constructor
  HF(const std::shared_ptr<const Input> &in_, const bool verbose_);
  explicit HF(const std::shared_ptr<const Input> &in_)
      : HF(in_, true) {}
  // Destructor
  virtual ~HF() = default;
  // Compute the scheme
  int compute();
  // Getters
  const Vector2D &getIdr() const { return idr; }
  const Vector2D &getLfc() const { return lfc; }
  const std::vector<double> &getSsf() const { return ssf; }
  const std::vector<double> &getWvg() const { return wvg; }
  std::vector<double> getRdf(const std::vector<double> &r) const;
  std::vector<double> getSdr() const;
  double getUInt() const;

protected:

  // Input data
  const std::shared_ptr<const Input> inPtr;
  // Integrator
  const std::shared_ptr<Integrator1D> itg;
  // Wave vector grid
  std::vector<double> wvg;
  // Ideal density response
  Vector2D idr;
  // Static local field correction
  Vector2D lfc;
  // Static structure factor
  std::vector<double> ssf;
  // Chemical potential
  double mu;
  // Access input pointer
  const Input &in() const { return *inPtr; }
  // Initialize basic properties
  virtual void init();
  // Calculations to compute the structural properties
  virtual void computeStructuralProperties();
  // Compute static structure factor
  virtual void computeSsf();
  virtual void computeSsfFinite();
  virtual void computeSsfGround();
  // Compute local field correction
  virtual void computeLfc();

private:

  // Construct wave vector grid
  void buildWaveVectorGrid();
  // Compute chemical potential
  void computeChemicalPotential();
  // Compute the ideal density response
  void computeIdr();
  void computeIdrFinite();
  void computeIdrGround();
};

namespace HFUtil {

  class Idr : public dimensionsUtil::DimensionsHandler {

  public:

    // Constructor remains unchanged
    Idr(const std::shared_ptr<const Input> in_,
        const double &x_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        std::shared_ptr<Integrator1D> itg_)
        : in(in_),
          x(x_),
          mu(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_),
          res(in_->getNMatsubara()) {}

    // Get result of integration
    std::vector<double> get();

  private:

    // Input parameters
    const std::shared_ptr<const Input> in;
    // Wave-vector
    const double x;
    // Chemical potential
    const double mu;
    // Integration limits for finite temperature calculations
    const double yMin;
    const double yMax;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Vector to store results
    std::vector<double> res;
    // Compute for 3D systems
    void compute3D() override;
    // Compute for 2D systems
    void compute2D() override;
    // Idr integrand for frequency = l and wave-vector x
    double integrand(const double &y, const int &l) const;
    // Idr integrand for frequency = 0 and wave-vector x
    double integrand(const double &y) const;
    // Idr integrand for 2D frequency = l and wave-vector x
    double integrand2D(const double &y, const int &l) const;
    // Idr integrand for 2D frequency = 0 and wave-vector x
    double integrand2D(const double &y) const;

  };

  class IdrGround {

  public:

    // Constructor
    IdrGround(const double &x_, const double &Omega_)
        : x(x_),
          Omega(Omega_) {}
    // Get
    double get() const;

  private:

    // Wave-vector
    const double x;
    // Frequency
    const double Omega;
  };

  class Ssf : public dimensionsUtil::DimensionsHandler {

  public:

    // Constructor for finite temperature calculations
    Ssf(const std::shared_ptr<const Input> in_,
        const double &x_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        std::shared_ptr<Integrator1D> itg_,
        const std::vector<double> &itgGrid_,
        std::shared_ptr<Integrator2D> itg2_,
        const Vector2D &idr_,
        const double &grid_val_)
        : in(in_),
          x(x_),
          mu(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_),
          itgGrid(itgGrid_),
          itg2(itg2_),
          idr(idr_),
          grid_val(grid_val_),
          res(x_) {}

    // Get at any temperature
    double get();

  private:
    // Input parameters
    const std::shared_ptr<const Input> in;
    // Wave-vector
    const double x;
    // Chemical potential
    const double mu;
    // Integration limits for finite temperature calculations
    const double yMin;
    const double yMax;
    void compute3D() override;
    // Compute for 2D systems
    void compute2D() override;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg2;
    const Vector2D idr;
    const double grid_val;
    double res;
    // Get integrand
    double integrand(const double &y) const;
    double integrand2DOut(const double &y) const;
    double integrand2DIn(const double &p) const;
    // Get at zero temperature
    double get0() const;
  };

  class SsfGround {

  public:

    // Constructor for zero temperature calculations
    explicit SsfGround(const double &x_)
        : x(x_) {}
    // Get result
    double get() const;

  private:

    // Wave-vector
    const double x;
  };

} // namespace HFUtil

#endif