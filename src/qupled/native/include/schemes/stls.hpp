#ifndef STLS_HPP
#define STLS_HPP

#include "schemes/input.hpp"
#include "schemes/rpa.hpp"
#include "util/mpi_util.hpp"
#include "util/numerics.hpp"
#include <cmath>
#include <vector>

/**
 * @brief Solver for the STLS (Singwi–Tosi–Land–Sjölander) dielectric scheme.
 *
 * Extends the RPA solver by computing the static local field correction (SLFC)
 * self-consistently via an iterative scheme. The SLFC accounts for short-range
 * correlations beyond the random phase approximation.
 */
class Stls : public Rpa {

public:

  /**
   * @brief Construct with explicit verbosity flag.
   * @param in_      Shared pointer to the STLS input parameters.
   * @param verbose_ If false, solver output is suppressed.
   */
  Stls(const std::shared_ptr<const StlsInput> &in_, const bool verbose_);

  /**
   * @brief Construct with verbosity enabled.
   * @param in_ Shared pointer to the STLS input parameters.
   */
  explicit Stls(const std::shared_ptr<const StlsInput> &in_)
      : Stls(in_, true) {}

  /** @brief Virtual destructor. */
  ~Stls() override = default;

  /**
   * @brief Return the convergence error from the last iteration.
   * @return RMS difference between successive SSF iterates.
   */
  double getError() const { return computeError(); }

protected:

  /** @brief SSF from the previous iteration (used for mixing). */
  std::vector<double> ssfOld;

  /** @brief Compute the SSF and SLFC self-consistently. */
  void computeStructuralProperties() override;

  /** @brief Compute the static structure factor using the current SLFC. */
  void computeSsf() override;

  /** @brief Compute the static local field correction from the current SSF. */
  void computeLfc() override;

  /** @brief Set the initial guess for the iterative procedure. */
  virtual void initialGuess();

  /**
   * @brief Attempt to load the initial guess from the input parameters.
   * @return True if a valid guess was found in the input, false otherwise.
   */
  virtual bool initialGuessFromInput();

  /**
   * @brief Compute the convergence error between successive iterates.
   * @return RMS difference between @p ssf and @p ssfOld.
   */
  virtual double computeError() const;

  /** @brief Mix the new and old SSF to advance the iteration. */
  virtual void updateSolution();

private:

  /** @brief Access the input as an @p StlsInput reference. */
  const StlsInput &in() const;
};

/** @brief Internal helpers for the STLS static local field correction. */
namespace StlsUtil {

  /**
   * @brief Perform a checked dynamic pointer cast between shared_ptr types.
   *
   * Throws an MPI error if the cast result is null (i.e., the runtime type
   * does not match @p TOut).
   *
   * @tparam T   Resulting shared_ptr type.
   * @param  ptr Pointer to validate.
   * @return @p ptr unchanged.
   */
  template <typename T>
  T check_dynamic_cast_result(T ptr) {
    if (!ptr) { MPIUtil::throwError("Unable to perform dynamic cast"); }
    return ptr;
  }

  /**
   * @brief Downcast a @c shared_ptr<const TIn> to @c shared_ptr<const TOut>.
   *
   * @tparam TIn  Source type.
   * @tparam TOut Target type (must be derived from or base of @p TIn).
   * @param  in   Source pointer.
   * @return Downcast pointer; throws on failure.
   */
  template <typename TIn, typename TOut>
  std::shared_ptr<const TOut>
  dynamic_pointer_cast(const std::shared_ptr<const TIn> &in) {
    return check_dynamic_cast_result(std::dynamic_pointer_cast<const TOut>(in));
  }

  /**
   * @brief Downcast a @c shared_ptr<TIn> to @c shared_ptr<TOut>.
   *
   * @tparam TIn  Source type.
   * @tparam TOut Target type (must be derived from or base of @p TIn).
   * @param  in   Source pointer.
   * @return Downcast pointer; throws on failure.
   */
  template <typename TIn, typename TOut>
  std::shared_ptr<TOut> dynamic_pointer_cast(const std::shared_ptr<TIn> &in) {
    return check_dynamic_cast_result(std::dynamic_pointer_cast<TOut>(in));
  }

  /**
   * @brief Base class holding shared state for SLFC helpers.
   */
  class SlfcBase {

  protected:

    /**
     * @brief Construct with the quantities needed for SLFC evaluation.
     * @param x_    Wave-vector value.
     * @param yMin_ Lower integration limit.
     * @param yMax_ Upper integration limit.
     * @param ssfi_ Shared pointer to an SSF interpolator.
     */
    SlfcBase(const double &x_,
             const double &yMin_,
             const double &yMax_,
             std::shared_ptr<Interpolator1D> ssfi_)
        : x(x_),
          yMin(yMin_),
          yMax(yMax_),
          ssfi(ssfi_) {}

    /** @brief Wave-vector value. */
    const double x;
    /** @brief Lower integration limit. */
    const double yMin;
    /** @brief Upper integration limit. */
    const double yMax;
    /** @brief Interpolator for the static structure factor. */
    const std::shared_ptr<Interpolator1D> ssfi;

    /**
     * @brief Evaluate the interpolated SSF at auxiliary momentum @p y.
     * @param y Auxiliary momentum value.
     * @return Interpolated SSF value.
     */
    double ssf(const double &y) const;
  };

  /**
   * @brief Computes the STLS static local field correction.
   *
   * Evaluates the SLFC integral at a single wave-vector @p x_ by integrating
   * over the auxiliary momentum with the SSF interpolator.
   */
  class Slfc : public SlfcBase, dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a SLFC calculation.
     * @param x_    Wave-vector value.
     * @param yMin_ Lower integration limit.
     * @param yMax_ Upper integration limit.
     * @param ssfi_ Shared pointer to an SSF interpolator.
     * @param itg_  Shared pointer to a 1D integrator.
     * @param in_   Shared pointer to the input parameters.
     */
    Slfc(const double &x_,
         const double &yMin_,
         const double &yMax_,
         std::shared_ptr<Interpolator1D> ssfi_,
         std::shared_ptr<Integrator1D> itg_,
         const std::shared_ptr<const Input> in_)
        : SlfcBase(x_, yMin_, yMax_, ssfi_),
          itg(itg_),
          in(in_),
          res(numUtil::NaN) {}

    /**
     * @brief Compute and return the SLFC at the current wave-vector.
     * @return SLFC value.
     */
    double get();

  private:

    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /**
     * @brief 3D integrand over auxiliary momentum @p y.
     * @param y Auxiliary momentum variable.
     */
    double integrand(const double &y) const;
    /**
     * @brief 2D integrand over auxiliary momentum @p y.
     * @param y Auxiliary momentum variable.
     */
    double integrand2D(const double &y) const;
    /** @brief Stores the result of the integration. */
    double res;
    void compute2D() override;
    void compute3D() override;
  };

} // namespace StlsUtil

#endif
