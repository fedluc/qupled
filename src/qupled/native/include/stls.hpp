#ifndef STLS_HPP
#define STLS_HPP

#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "rpa.hpp"
#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Solver for the STLS scheme
// -----------------------------------------------------------------

class Stls : public Rpa {

public:

  // Constructors
  Stls(const std::shared_ptr<const StlsInput> &in_, const bool verbose_);
  explicit Stls(const std::shared_ptr<const StlsInput> &in_)
      : Stls(in_, true) {}
  // Destructor
  ~Stls() override = default;
  // Getters
  double getError() const { return computeError(); }

protected:

  // Static local field correction to use during the iterations
  std::vector<double> ssfOld;
  // Compute structural properties
  void computeStructuralProperties() override;
  // Compute static structure factor
  void computeSsf() override;
  // Compute static local field correction
  void computeLfc() override;
  // Iterations to solve the stls scheme
  void initialGuess();
  virtual bool initialGuessFromInput();
  double computeError() const;
  virtual void updateSolution();

private:

  // Input parameters
  const StlsInput &in() const;
};

namespace StlsUtil {

  // --------------------------------------------------------
  // Method to dyanmic cast between different shared pointers
  // --------------------------------------------------------

  template <typename T>
  T check_dynamic_cast_result(T ptr) {
    if (!ptr) { MPIUtil::throwError("Unable to perform dynamic cast"); }
    return ptr;
  }

  template <typename TIn, typename TOut>
  std::shared_ptr<const TOut>
  dynamic_pointer_cast(const std::shared_ptr<const TIn> &in) {
    return check_dynamic_cast_result(std::dynamic_pointer_cast<const TOut>(in));
  }

  template <typename TIn, typename TOut>
  std::shared_ptr<TOut> dynamic_pointer_cast(const std::shared_ptr<TIn> &in) {
    return check_dynamic_cast_result(std::dynamic_pointer_cast<TOut>(in));
  }

  // -----------------------------------------------------------------
  // Classes for the static local field correction
  // -----------------------------------------------------------------

  class SlfcBase {

  protected:

    // Constructor
    SlfcBase(const double &x_,
             const double &yMin_,
             const double &yMax_,
             std::shared_ptr<Interpolator1D> ssfi_)
        : x(x_),
          yMin(yMin_),
          yMax(yMax_),
          ssfi(ssfi_) {}
    // Wave-vector
    const double x;
    // Integration limits
    const double yMin;
    const double yMax;
    // Static structure factor interpolator
    const std::shared_ptr<Interpolator1D> ssfi;
    // Compute static structure factor
    double ssf(const double &y) const;
  };

  class Slfc : public SlfcBase, dimensionsUtil::DimensionsHandler {

  public:

    // Constructor
    Slfc(const double &x_,
         const double &yMin_,
         const double &yMax_,
         std::shared_ptr<Interpolator1D> ssfi_,
         std::shared_ptr<Integrator1D> itg_,
         const std::shared_ptr<const Input> in_)
        : SlfcBase(x_, yMin_, yMax_, ssfi_),
          itg(itg_), in(in_), res(numUtil::NaN) {}
    // Get result of integration
    double get();

  private:

    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    const std::shared_ptr<const Input> in;
    // Integrand
    double integrand(const double &y) const;
    double integrand2D(const double &y) const;
    // result of integration
    double res;
    void compute2D() override;
    void compute3D() override;

  };

} // namespace StlsUtil

#endif
