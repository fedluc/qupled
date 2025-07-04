#ifndef HF2D_HPP
#define HF2D_HPP

#include "hf.hpp"
#include "numerics.hpp"
#include <vector>

// -----------------------------------------------------------------
// Solver for the 2D Hartree-Fock scheme
// -----------------------------------------------------------------

class HF2D : public HF {

public:

  // Constructor
  HF2D(const std::shared_ptr<const Input> &in_, const bool verbose_);
  explicit HF2D(const std::shared_ptr<const Input> &in_)
      : HF2D(in_, true) {}
  // Destructor
  virtual ~HF2D() = default;
  // Compute the scheme
  int compute2D();
  // Getters
  const Vector2D &getIdr2D() const { return idr2D; }
  const Vector2D &getLfc2D() const { return lfc2D; }
  const std::vector<double> &getSsf2D() const { return ssf2D; }
  const std::vector<double> &getWvg() const { return wvg; }
  std::vector<double> getRdf2D(const std::vector<double> &r) const;
  std::vector<double> getSdr2D() const;
  double getUInt2D() const;

protected:

  // Input data
  const std::shared_ptr<const Input> inPtr;
  // Integrator
  const std::shared_ptr<Integrator1D> itg;
  const std::shared_ptr<Integrator2D> itg2;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  // Wave vector grid
  std::vector<double> wvg;
  // Ideal density response
  Vector2D idr2D;
  // Static local field correction
  Vector2D lfc2D;
  // Static structure factor
  std::vector<double> ssf2D;
  // Chemical potential
  double mu2D;
  // Access input pointer
  const Input &in() const { return *inPtr; }
  // Initialize basic properties
  virtual void init2D();
  // Calculations to compute the structural properties
  virtual void computeStructuralProperties2D();
  // Compute static structure factor
  virtual void computeSsf2D();
  virtual void computeSsfFinite2D();
  virtual void computeSsfGround2D();
  // Compute local field correction
  virtual void computeLfc2D();

private:

  // Construct wave vector grid
  void buildWaveVectorGrid();
  // Compute chemical potential
  void computeChemicalPotential2D();
  // Compute the ideal density response
  void computeIdr2D();
  void computeIdrFinite2D();
  void computeIdrGround2D();
};

namespace HFUtil2D {

  class Idr2D {

  public:

    // Constructor
    Idr2D(const int nl_,
        const double &x_,
        const double &Theta_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        std::shared_ptr<Integrator1D> itg_)
        : nl(nl_),
          x(x_),
          Theta(Theta_),
          mu2D(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_) {}
    // Get result of integration
    std::vector<double> get2D() const;

  private:

    // Number of matsubara frequency
    const int nl;
    // Wave-vector
    const double x;
    // Degeneracy parameter
    const double Theta;
    // Chemical potential
    const double mu2D;
    // Integration limits for finite temperature calculations
    const double yMin;
    const double yMax;
    // Idr integrand for frequency = l and wave-vector x
    double integrand2D(const double &y, const int &l) const;
    // Idr integrand for frequency = 0 and wave-vector x
    double integrand2D(const double &y) const;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
  };

  class IdrGround2D {

  public:

    // Constructor
    IdrGround2D(const double &x_, const double &Omega_)
        : x(x_),
          Omega(Omega_) {}
    // Get
    double get2D() const;

  private:

    // Wave-vector
    const double x;
    // Frequency
    const double Omega;
  };

  class Ssf2D {

  public:

    // Constructor for finite temperature calculations
    Ssf2D(const double &x_,
        const double &Theta_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        const std::vector<double> &itgGrid_,
        std::shared_ptr<Integrator2D> itg2_)
        : x(x_),
          Theta(Theta_),
          mu2D(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itgGrid(itgGrid_),
          itg2(itg2_) {}
    // Get at any temperature
    double get2D() const;

  private:

    // Wave-vector
    const double x;
    // Degeneracy parameter
    const double Theta;
    // Chemical potential
    const double mu2D;
    // Integration limits for finite temperature calculations
    const double yMin;
    const double yMax;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg2;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
    // Get integrand
    double integrandOut(const double &y) const;
    double integrandIn(const double &p) const;
    // Get at zero temperature
    double get0() const;
  };

  class SsfGround2D {

  public:

    // Constructor for zero temperature calculations
    explicit SsfGround2D(const double &x_)
        : x(x_) {}
    // Get result
    double get() const;

  private:

    // Wave-vector
    const double x;
  };

} // namespace HFUtil

#endif