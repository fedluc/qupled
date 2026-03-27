#ifndef QSTLS_HPP
#define QSTLS_HPP

#include "schemes/input.hpp"
#include "schemes/stls.hpp"
#include "util/numerics.hpp"
#include "util/vector2D.hpp"
#include "util/vector3D.hpp"
#include <map>

/**
 * @brief Solver for the qSTLS (quantum STLS) dielectric scheme.
 *
 * Extends the STLS solver by replacing the classical SLFC with a quantum
 * local field correction derived from the auxiliary density response (ADR).
 * The ADR is split into a "fixed" frequency-independent component (stored in
 * an SQLite database) and a frequency-dependent part that is iterated to
 * self-consistency together with the SSF.
 */
class Qstls : public Stls {

public:

  /**
   * @brief Construct with explicit verbosity flag.
   * @param in_      Shared pointer to the qSTLS input parameters.
   * @param verbose_ If false, solver output is suppressed.
   */
  Qstls(const std::shared_ptr<const QstlsInput> &in_, const bool verbose_);

  /**
   * @brief Construct with verbosity enabled.
   * @param in_ Shared pointer to the qSTLS input parameters.
   */
  explicit Qstls(const std::shared_ptr<const QstlsInput> &in_)
      : Qstls(in_, true) {}

  /**
   * @brief Return the fixed component of the auxiliary density response.
   * @return 3D array indexed by (wave-vector, wave-vector, Matsubara
   * frequency).
   */
  const Vector3D &getAdrFixed() const { return adrFixed; }

protected:

  /** @brief Fixed (frequency-independent) component of the ADR. */
  Vector3D adrFixed;

  /** @brief Path to the SQLite database used to persist @p adrFixed. */
  std::string adrFixedDatabaseName;

  /** @brief Initialize wave-vector grid, chemical potential, and ADR fixed
   * component. */
  void init() override;

  /** @brief Compute the quantum local field correction from the ADR. */
  void computeLfc() override;

  /**
   * @brief Read the fixed ADR component from the database.
   * @param res   Output array to populate.
   * @param name  Database file path.
   * @param runId Run identifier in the database.
   */
  void readAdrFixed(Vector3D &res, const std::string &name, int runId) const;

  /**
   * @brief Write the fixed ADR component to the database.
   * @param res  Array to persist.
   * @param name Database file path.
   */
  void writeAdrFixed(const Vector3D &res, const std::string &name) const;

private:

  /** @brief Access the input as a @p QstlsInput reference. */
  const QstlsInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, QstlsInput>(inPtr);
  }

  /** @brief Compute and cache the fixed component of the ADR. */
  void computeAdrFixed();

  /** @brief Compute the zero-temperature static structure factor. */
  void computeSsfGround() override;
};

/** @brief Internal helpers for the qSTLS auxiliary density response. */
namespace QstlsUtil {

  /// @cond INTERNAL
  constexpr const char *SQL_TABLE_NAME = "fixed";
  constexpr const char *SQL_CREATE_TABLE = R"(
      CREATE TABLE IF NOT EXISTS {} (
          run_id INTEGER NOT NULL,
          name TEXT NOT NULL,
          value TEXT NOT NULL,
          PRIMARY KEY (run_id, name),
          FOREIGN KEY(run_id) REFERENCES {}(id) ON DELETE CASCADE
      );
    )";
  constexpr const char *SQL_INSERT =
      "INSERT OR REPLACE INTO {} (run_id, name, value) VALUES (?, ?, ?);";
  constexpr const char *SQL_SELECT =
      "SELECT value FROM {} WHERE run_id = ? AND name = ?;";
  constexpr const char *SQL_SELECT_RUN_ID =
      "SELECT value FROM {} WHERE run_id = ?";
  constexpr const char *SQL_SELECT_TABLE =
      "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?";
  /// @endcond

  /**
   * @brief Delete blob data associated with a specific run from disk.
   * @param dbName  Path to the SQLite database file.
   * @param runId   Run identifier whose blobs should be deleted.
   */
  void deleteBlobDataOnDisk(const std::string &dbName, int runId);

  /**
   * @brief Base class holding shared state for ADR helper classes.
   */
  class AdrBase {

  public:

    /**
     * @brief Construct with thermodynamic and grid parameters.
     * @param Theta_ Degeneracy parameter.
     * @param yMin_  Lower integration limit.
     * @param yMax_  Upper integration limit.
     * @param x_     Wave-vector value.
     * @param ssfi_  Shared pointer to an SSF interpolator.
     */
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

    /** @brief Degeneracy parameter. */
    const double Theta;
    /** @brief Lower integration limit. */
    const double yMin;
    /** @brief Upper integration limit. */
    const double yMax;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Interpolator for the static structure factor. */
    const std::shared_ptr<Interpolator1D> ssfi;
    /** @brief Integrand scaling constant. */
    const double isc;
    /** @brief Integrand scaling constant for the zeroth Matsubara frequency. */
    const double isc0;

    /**
     * @brief Evaluate the interpolated SSF at auxiliary momentum @p y.
     * @param y Auxiliary momentum.
     * @return Interpolated SSF value.
     */
    double ssf(const double &y) const;
  };

  /**
   * @brief Base class holding shared state for fixed-ADR helper classes.
   */
  class AdrFixedBase {

  public:

    /**
     * @brief Construct with thermodynamic and grid parameters.
     * @param Theta_ Degeneracy parameter.
     * @param qMin_  Lower integration limit over the auxiliary momentum.
     * @param qMax_  Upper integration limit over the auxiliary momentum.
     * @param x_     Wave-vector value.
     * @param mu_    Chemical potential.
     */
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

    /** @brief Degeneracy parameter. */
    const double Theta;
    /** @brief Lower integration limit. */
    const double qMin;
    /** @brief Upper integration limit. */
    const double qMax;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Chemical potential. */
    const double mu;
  };

  /**
   * @brief Computes the finite-temperature auxiliary density response.
   *
   * Evaluates the ADR at wave-vector @p x_ and all Matsubara frequencies by
   * integrating over the auxiliary momentum using the fixed ADR component and
   * the current SSF.
   */
  class Adr : public AdrBase {

  public:

    /**
     * @brief Construct for a finite-temperature ADR calculation.
     * @param Theta_ Degeneracy parameter.
     * @param yMin_  Lower integration limit.
     * @param yMax_  Upper integration limit.
     * @param x_     Wave-vector value.
     * @param ssfi_  Shared pointer to an SSF interpolator.
     * @param itg_   Shared pointer to a 1D integrator.
     */
    Adr(const double &Theta_,
        const double &yMin_,
        const double &yMax_,
        const double &x_,
        std::shared_ptr<Interpolator1D> ssfi_,
        std::shared_ptr<Integrator1D> itg_)
        : AdrBase(Theta_, yMin_, yMax_, x_, ssfi_),
          itg(itg_) {}

    /**
     * @brief Compute and store the ADR for all Matsubara frequencies.
     * @param wvg   Wave-vector grid.
     * @param fixed Fixed ADR component (wave-vectors × wave-vectors ×
     * frequencies).
     * @param res   Output array (wave-vectors × frequencies) to accumulate
     * into.
     */
    void
    get(const std::vector<double> &wvg, const Vector3D &fixed, Vector2D &res);

  private:

    /**
     * @brief Evaluate the fixed-ADR contribution at auxiliary momentum @p y.
     * @param y Auxiliary momentum.
     * @return Fixed-component value at @p y.
     */
    double fix(const double &y) const;

    /**
     * @brief Integrand for the ADR frequency integral.
     * @param y Auxiliary momentum.
     * @return Integrand value at @p y.
     */
    double integrand(const double &y) const;

    /** @brief Interpolator for the fixed ADR component. */
    Interpolator1D fixi;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
  };

  /**
   * @brief Computes the fixed (frequency-independent) component of the ADR.
   *
   * This component depends only on the Fermi-Dirac occupation numbers and does
   * not change during the self-consistency iterations.
   */
  class AdrFixed : public AdrFixedBase {

  public:

    /**
     * @brief Construct for a finite-temperature fixed-ADR calculation.
     * @param Theta_    Degeneracy parameter.
     * @param qMin_     Lower integration limit.
     * @param qMax_     Upper integration limit.
     * @param x_        Wave-vector value.
     * @param mu_       Chemical potential.
     * @param itgGrid_  Grid for the 2D integration.
     * @param itg_      Shared pointer to a 2D integrator.
     */
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

    /**
     * @brief Compute and store the fixed ADR component for all wave-vector
     * pairs.
     * @param wvg Wave-vector grid.
     * @param res Output 3D array (wave-vectors × wave-vectors × frequencies).
     */
    void get(const std::vector<double> &wvg, Vector3D &res) const;

  private:

    /**
     * @brief First integrand over auxiliary momentum @p q and frequency @p l.
     * @param q Auxiliary momentum.
     * @param l Matsubara frequency value.
     */
    double integrand1(const double &q, const double &l) const;

    /**
     * @brief Second integrand over variables @p t, @p y, and @p l.
     * @param t Intermediate momentum variable.
     * @param y Auxiliary momentum variable.
     * @param l Matsubara frequency value.
     */
    double integrand2(const double &t, const double &y, const double &l) const;

    /** @brief 2D numerical integrator. */
    const std::shared_ptr<Integrator2D> itg;
    /** @brief Grid for 2D integration. */
    const std::vector<double> &itgGrid;
  };

  /**
   * @brief Computes the auxiliary density response at zero temperature.
   *
   * Evaluates the ground-state ADR by integrating over the real frequency
   * and the auxiliary momentum.
   */
  class AdrGround : public AdrBase {

  public:

    /**
     * @brief Construct for a zero-temperature ADR calculation.
     * @param x_     Wave-vector value.
     * @param Omega_ Real frequency value.
     * @param ssfi_  Shared pointer to an SSF interpolator.
     * @param yMax_  Upper integration limit.
     * @param itg_   Shared pointer to a 2D integrator.
     */
    AdrGround(const double &x_,
              const double &Omega_,
              std::shared_ptr<Interpolator1D> ssfi_,
              const double &yMax_,
              std::shared_ptr<Integrator2D> itg_)
        : AdrBase(0.0, 0.0, yMax_, x_, ssfi_),
          Omega(Omega_),
          itg(itg_) {}

    /**
     * @brief Compute and return the ground-state ADR.
     * @return ADR value at (@p x_, @p Omega_).
     */
    double get();

  private:

    /** @brief Real frequency. */
    const double Omega;
    /** @brief 2D numerical integrator. */
    const std::shared_ptr<Integrator2D> itg;
    /**
     * @brief Outer integrand over auxiliary momentum @p y.
     * @param y Auxiliary momentum.
     */
    double integrand1(const double &y) const;
    /**
     * @brief Inner integrand over intermediate variable @p t.
     * @param t Intermediate integration variable.
     */
    double integrand2(const double &t) const;
  };

  /**
   * @brief Computes the zero-temperature qSTLS static structure factor.
   *
   * Extends the RPA ground-state SSF to include the quantum local field
   * correction derived from the ADR.
   */
  class SsfGround : public RpaUtil::SsfGround {

  public:

    /**
     * @brief Construct for a zero-temperature qSTLS SSF calculation.
     * @param x_     Wave-vector value.
     * @param ssfHF_ Hartree-Fock SSF at @p x_.
     * @param xMax_  Upper integration limit for the wave-vector integral.
     * @param ssfi_  Shared pointer to an SSF interpolator.
     * @param itg_   Shared pointer to a 1D integrator.
     * @param in_    Shared pointer to the input parameters.
     */
    SsfGround(const double &x_,
              const double &ssfHF_,
              const double &xMax_,
              std::shared_ptr<Interpolator1D> ssfi_,
              std::shared_ptr<Integrator1D> itg_,
              const std::shared_ptr<const Input> in_)
        : RpaUtil::SsfGround(x_, ssfHF_, std::span<const double>(), itg_, in_),
          xMax(xMax_),
          ssfi(ssfi_) {}

    /**
     * @brief Compute and return the zero-temperature qSTLS SSF.
     * @return SSF value at @p x_.
     */
    double get();

  private:

    /** @brief Upper integration limit for the wave-vector integral. */
    const double xMax;
    /** @brief SSF interpolator. */
    const std::shared_ptr<Interpolator1D> ssfi;
    /**
     * @brief Integrand for the zero-temperature frequency integral.
     * @param Omega Real frequency value.
     * @return Integrand value at @p Omega.
     */
    double integrand(const double &Omega) const;
  };

} // namespace QstlsUtil

#endif
