#ifndef RECSTLS_HPP
#define RECSTLS_HPP

#include "numerics.hpp"
#include "rpa.hpp"

class RecStls : public Rpa {

public:

  RecStls(const std::shared_ptr<const RecStlsInput> &in_, const bool verbose_);
  explicit RecStls(const std::shared_ptr<const RecStlsInput> &in_)
      : RecStls(in_, true) {}
  ~RecStls() override = default;

  const std::vector<double> &getSsfInput() const { return ssfInput; }

protected:

  void computeStructuralProperties() override;
  void computeLfc() override;

private:

  const RecStlsInput &in() const;
  void computeInputSsf();

  std::vector<double> ssfInput;
  std::shared_ptr<Interpolator1D> rdfi;
};

namespace RecStlsUtil {

  class InputSsf {

  public:

    InputSsf(const double &x_,
             const double &uMin_,
             const double &uMax_,
             std::shared_ptr<Interpolator1D> rdfi_,
             std::shared_ptr<Integrator1D> itg_)
        : x(x_),
          uMin(uMin_),
          uMax(uMax_),
          rdfi(rdfi_),
          itg(itg_) {}

    double get() const;

  private:

    const double x;
    const double uMin;
    const double uMax;
    const std::shared_ptr<Interpolator1D> rdfi;
    const std::shared_ptr<Integrator1D> itg;

    double integrand(const double &u) const;
  };

} // namespace RecStlsUtil

#endif
