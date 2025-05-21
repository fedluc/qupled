#ifndef QSTLS_HPP
#define QSTLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "stls.hpp"
#include "vector2D.hpp"
#include "vector3D.hpp"
#include <map>

// -----------------------------------------------------------------
// Solver for the qSTLS-based schemes
// -----------------------------------------------------------------

class Qstls : public Stls {

public:

  // Constructor
  Qstls(const std::shared_ptr<const QstlsInput> &in_, const bool verbose_);
  Qstls(const std::shared_ptr<const QstlsInput> &in_)
      : Qstls(in_, true) {}
  Qstls(const QstlsInput &in_, const bool verbose_)
      : Qstls(std::make_shared<const QstlsInput>(in_), verbose_) {}
  explicit Qstls(const QstlsInput &in_)
      : Qstls(in_, true) {}

  // Getters
  const Vector3D &getAdrFixed() const { return adrFixed; }

protected:

  // Auxiliary density response
  Vector3D adrFixed;
  std::string adrFixedDatabaseName;
  // Initialize basic properties
  void init() override;
  // Compute auxiliary density response
  void computeLfc() override;
  // Read and write auxiliary density response to database
  void readAdrFixed(Vector3D &res, const std::string &name, int runId) const;
  void writeAdrFixed(const Vector3D &res, const std::string &name) const;

private:

  // Input parameters
  const QstlsInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, QstlsInput>(inPtr);
  }
  // Compute auxiliary density response
  void computeAdrFixed();
  // Compute static structure factor at finite temperature
  void computeSsfGround() override;
};

namespace QstlsUtil {

  // -----------------------------------------------------------------
  // SQL queries
  // -----------------------------------------------------------------

  constexpr const char *SQL_TABLE_NAME = "fixed";

  constexpr const char *SQL_CREATE_TABLE = R"(
      CREATE TABLE IF NOT EXISTS {} (
          run_id INTEGER NOT NULL,
          name TEXT NOT NULL,
          adr BLOB NOT NULL,
          PRIMARY KEY (run_id, name),
          FOREIGN KEY(run_id) REFERENCES {}(id) ON DELETE CASCADE
      );
    )";

  constexpr const char *SQL_INSERT =
      "INSERT OR REPLACE INTO {} (run_id, name, adr) VALUES (?, ?, ?);";

  constexpr const char *SQL_SELECT =
      "SELECT adr FROM {} WHERE run_id = ? AND name = ?;";
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
            std::shared_ptr<Interpolator1D> ssfi_)
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
    const std::shared_ptr<Interpolator1D> ssfi;
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
        std::shared_ptr<Interpolator1D> ssfi_,
        std::shared_ptr<Integrator1D> itg_)
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
    const std::shared_ptr<Integrator1D> itg;
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
             std::shared_ptr<Integrator2D> itg_)
        : AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
          itg(itg_),
          itgGrid(itgGrid_) {}

    // Get integration result
    void get(const std::vector<double> &wvg, Vector3D &res) const;

  private:

    // Integrands
    double integrand1(const double &q, const double &l) const;
    double integrand2(const double &t, const double &y, const double &l) const;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
  };

  class AdrGround : public AdrBase {

  public:

    // Constructor for zero temperature calculations
    AdrGround(const double &x_,
              const double &Omega_,
              std::shared_ptr<Interpolator1D> ssfi_,
              const double &yMax_,
              std::shared_ptr<Integrator2D> itg_)
        : AdrBase(0.0, 0.0, yMax_, x_, ssfi_),
          Omega(Omega_),
          itg(itg_) {}
    // Get
    double get();

  private:

    // Frequency
    const double Omega;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg;
    // Integrands
    double integrand1(const double &y) const;
    double integrand2(const double &t) const;
  };

  // -----------------------------------------------------------------
  // Class for the static structure factor
  // -----------------------------------------------------------------

  class SsfGround : public RpaUtil::SsfGround {

  public:

    // Constructor for zero temperature calculations
    SsfGround(const double &x_,
              const double &rs_,
              const double &ssfHF_,
              const double &xMax_,
              const double &OmegaMax_,
              std::shared_ptr<Interpolator1D> ssfi_,
              std::shared_ptr<Integrator1D> itg_)
        : RpaUtil::SsfGround(
              x_, rs_, ssfHF_, std::span<const double>(), OmegaMax_, itg_),
          xMax(xMax_),
          ssfi(ssfi_) {}
    // Get result of integration
    double get();

  private:

    // Integration limit for the wave-vector integral
    const double xMax;
    // Interpolator
    const std::shared_ptr<Interpolator1D> ssfi;
    // Integrand for zero temperature calculations
    double integrand(const double &Omega) const;
  };

} // namespace QstlsUtil

#endif
