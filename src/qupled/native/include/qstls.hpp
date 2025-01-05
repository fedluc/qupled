#ifndef QSTLS_HPP
#define QSTLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "stls.hpp"
#include "vector2D.hpp"
#include "vector3D.hpp"
#include <fmt/core.h>
#include <map>

// -----------------------------------------------------------------
// Solver for the qSTLS-based schemes
// -----------------------------------------------------------------

class Qstls : public Stls {

public:

  // Constructor
  Qstls(const QstlsInput &in_, const bool verbose_, const bool writeFiles_);
  explicit Qstls(const QstlsInput &in_)
      : Qstls(in_, true, true) {}
  // Compute qstls scheme
  int compute();
  // Getters
  const QstlsInput &getInput() const { return in; }
  double getError() const { return computeError(); }
  const Vector2D &getAdr() const { return adr; }
  const Vector3D &getAdrFixed() const { return adrFixed; }

protected:

  // Input data
  const QstlsInput in;
  // Auxiliary density response
  Vector2D adr;
  Vector2D adrOld;
  Vector3D adrFixed;
  std::string adrFixedFileName =
      fmt::format("adr_fixed_theta{:.3f}_matsubara{:}_{}.bin",
                  in.getDegeneracy(),
                  in.getNMatsubara(),
                  in.getTheory());
  std::map<int, std::pair<std::string, bool>> adrFixedIetFileInfo;
  // Static structure factor (for iterations)
  std::vector<double> ssfNew;
  std::vector<double> ssfOld;
  // Initialize basic properties
  void init();
  // Compute auxiliary density response
  void computeAdr();
  void computeAdrFixed();
  void writeAdrFixedFile(const Vector3D &res,
                         const std::string &fileName) const;
  void readAdrFixedFile(Vector3D &res,
                        const std::string &fileName,
                        const bool iet) const;
  bool checkAdrFixed(const std::vector<double> &wvg_,
                     const double Theta_,
                     const int nl_) const;
  void computeAdrIet();
  void computeAdrFixedIet();
  void getAdrFixedIetFileInfo();
  // Compute static structure factor at finite temperature
  void computeSsf();
  void computeSsfFinite();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  bool initialGuessFromRecovery();
  bool initialGuessFromInput();
  bool initialGuessSsf(const std::vector<double> &wvg_,
                       const std::vector<double> &adr_);
  bool initialGuessAdr(const std::vector<double> &wvg_, const Vector2D &adr_);
  bool initialGuessAdrFixed(const std::vector<double> &wvg_,
                            const double &Theta,
                            const int &nl_,
                            const Vector3D &adrFixed_);
  double computeError() const;
  void updateSolution();
  // Recovery files
  void writeRecovery();
  void readRecovery(std::vector<double> &wvg_,
                    std::vector<double> &ssf_,
                    Vector2D &adr_,
                    Vector3D &adrFixed_,
                    double &Theta,
                    int &nl) const;
};

// -----------------------------------------------------------------
// Class for the static structure factor
// -----------------------------------------------------------------

class Qssf : public Ssf {

public:

  // Constructor for quantum schemes
  Qssf(const double &x_,
       const double &Theta_,
       const double &rs_,
       const double &ssfHF_,
       const int &nl_,
       const double *idr_,
       const double *adr_,
       const double &bf_)
      : Ssf(x_, Theta_, rs_, ssfHF_, 0, nl_, idr_),
        adr(adr_),
        bf(bf_) {}
  // Get static structure factor
  double get() const;

private:

  // Auxiliary density response
  const double *adr;
  // Bridge function
  const double bf;
};
// -----------------------------------------------------------------
// Classes for the auxiliary density response
// -----------------------------------------------------------------

class AdrBase {

public:

  // Constructor
  AdrBase(const double &Theta_,
          const double &yMin_,
          const double &yMax_,
          const double &x_,
          const Interpolator1D &ssfi_)
      : Theta(Theta_),
        yMin(yMin_),
        yMax(yMax_),
        x(x_),
        ssfi(ssfi_),
        isc(-3.0 / 8.0),
        isc0(isc * 2.0 / Theta) {}

protected:

  // Degeneracy parameter
  const double Theta;
  // Integration limits
  const double yMin;
  const double yMax;
  // Wave-vector
  const double x;
  // Interpolator for the static structure factor
  const Interpolator1D &ssfi;
  // Integrand scaling constants
  const double isc;
  const double isc0;
  // Compute static structure factor
  double ssf(const double &y) const;
};

class AdrFixedBase {

public:

  // Constructor for finite temperature calculations
  AdrFixedBase(const double &Theta_,
               const double &qMin_,
               const double &qMax_,
               const double &x_,
               const double &mu_)
      : Theta(Theta_),
        qMin(qMin_),
        qMax(qMax_),
        x(x_),
        mu(mu_) {}

protected:

  // Degeneracy parameter
  const double Theta;
  // Integration limits
  const double qMin;
  const double qMax;
  // Wave-vector
  const double x;
  // Chemical potential
  const double mu;
};

class Adr : public AdrBase {

public:

  // Constructor for finite temperature calculations
  Adr(const double &Theta_,
      const double &yMin_,
      const double &yMax_,
      const double &x_,
      const Interpolator1D &ssfi_,
      Integrator1D &itg_)
      : AdrBase(Theta_, yMin_, yMax_, x_, ssfi_),
        itg(itg_) {}

  // Get result of integration
  void
  get(const std::vector<double> &wvg, const Vector3D &fixed, Vector2D &res);

private:

  // Compute fixed component
  double fix(const double &y) const;
  // integrand
  double integrand(const double &y) const;
  // Interpolator for the fixed component
  Interpolator1D fixi;
  // Integrator object
  Integrator1D &itg;
};

class AdrFixed : public AdrFixedBase {

public:

  // Constructor for finite temperature calculations
  AdrFixed(const double &Theta_,
           const double &qMin_,
           const double &qMax_,
           const double &x_,
           const double &mu_,
           const std::vector<double> &itgGrid_,
           Integrator2D &itg_)
      : AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
        itg(itg_),
        itgGrid(itgGrid_) {}

  // Get integration result
  void get(std::vector<double> &wvg, Vector3D &res) const;

private:

  // Integrands
  double integrand1(const double &q, const double &l) const;
  double integrand2(const double &t, const double &y, const double &l) const;
  // Integrator object
  Integrator2D &itg;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
};

// Class for the auxiliary density response calculation in the IET scheme
class AdrIet : public AdrBase {

public:

  // Constructor for finite temperature calculations
  AdrIet(const double &Theta_,
         const double &qMin_,
         const double &qMax_,
         const double &x_,
         const Interpolator1D &ssfi_,
         const std::vector<Interpolator1D> &dlfci_,
         const Interpolator1D &bfi_,
         const std::vector<double> &itgGrid_,
         Integrator2D &itg_)
      : AdrBase(Theta_, qMin_, qMax_, x_, ssfi_),
        itg(itg_),
        itgGrid(itgGrid_),
        dlfci(dlfci_),
        bfi(bfi_) {}

  // Get integration result
  void
  get(const std::vector<double> &wvg, const Vector3D &fixed, Vector2D &res);

private:

  // Integration limits
  const double &qMin = yMin;
  const double &qMax = yMax;
  // Integrands
  double integrand1(const double &q, const int &l) const;
  double integrand2(const double &y) const;
  // Integrator object
  Integrator2D &itg;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  // Interpolator for the dynamic local field correction
  const std::vector<Interpolator1D> &dlfci;
  // Interpolator for the bridge function contribution
  const Interpolator1D &bfi;
  // Interpolator for the fixed component
  Interpolator2D fixi;
  // Compute dynamic local field correction
  double dlfc(const double &y, const int &l) const;
  // Compute bridge function contribution
  double bf(const double &y) const;
  // Compute fixed component
  double fix(const double &x, const double &y) const;
};

class AdrFixedIet : public AdrFixedBase {

public:

  // Constructor for finite temperature calculations
  AdrFixedIet(const double &Theta_,
              const double &qMin_,
              const double &qMax_,
              const double &x_,
              const double &mu_,
              Integrator1D &itg_)
      : AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
        itg(itg_) {}

  // Get integration result
  void get(std::vector<double> &wvg, Vector3D &res) const;

private:

  // Integration limits
  const double &tMin = qMin;
  const double &tMax = qMax;
  // Integrands
  double integrand(const double &t,
                   const double &y,
                   const double &q,
                   const double &l) const;
  // Integrator object
  Integrator1D &itg;
};

// -----------------------------------------------------------------
// Ground state calculations
// -----------------------------------------------------------------

class AdrGround {

public:

  // Constructor for zero temperature calculations
  AdrGround(const double &Omega_,
            const double &x_,
            const Interpolator1D &ssfi_,
            const double &yMin_,
            const double &yMax_,
            Integrator1D &itg_)
      : Omega(Omega_),
        x(x_),
        ssfi(ssfi_),
        yMin(yMin_),
        yMax(yMax_),
        itg(itg_) {}
  // Get real part
  double real() const;
  // Get imaginary part
  double imag() const;

private:

  // Frequency
  const double Omega;
  // Wave-vector
  const double x;
  // Interpolator for the static structure factor
  const Interpolator1D &ssfi;
  // Integration limits for zero temperature calculations
  const double yMin;
  const double yMax;
  // Integrator object
  Integrator1D &itg;
  // Compute the static structure factor
  double ssf(const double &y) const;
  // Integrands
  double compute(const bool isReal) const;
  double integrand(const double &y, const bool isReal) const;

  // Auxiliary function
  class Gamma {
  public:

    Gamma(const bool isReal_)
        : isReal(isReal_) {}
    double get(const double &a, const double &b, const double &c) const;

  private:

    const bool isReal;
    double real(const double &a, const double &b, const double &c) const;
    double imag(const double &a, const double &b, const double &c) const;
    double Gamma1(const double &a, const double &b) const;
    double Gamma2(const double &a, const double &b, const double &c) const;
  };
};

class QSsfGround {

public:

  // Constructor for zero temperature calculations
  QSsfGround(const double &x_,
             const double &rs_,
             const double &xMax_,
             const double &ssfHF_,
             const Interpolator1D &ssfi_,
             Integrator1D &itg_)
      : x(x_),
        rs(rs_),
        ssfHF(ssfHF_),
        wMin(std::max(0.0, x * (x - 2.0))),
        wMax(x * (x + 2.0)),
        yMin(0.0),
        yMax(xMax_),
        itg(itg_),
        ssfi(ssfi_) {}
  // Get result of integration
  double get() const;

private:

  // Wave-vector
  const double x;
  // Coupling parameter
  const double rs;
  // HF static structure factor
  const double ssfHF;
  // Constant for unit conversion
  const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
  // Integration limits for zero temperature calculations
  const double wMin;
  const double wMax;
  const double yMin;
  const double yMax;
  // Integrator object
  Integrator1D &itg;
  // Interpolator
  const Interpolator1D &ssfi;
  // Interaction potential
  const double ip = 4.0 * lambda * rs / (M_PI * x * x);
  // Integrand for zero temperature calculations
  double integrand(const double &Omega) const;
};

#endif
