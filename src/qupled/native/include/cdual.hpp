#ifndef CDUAL_HPP
#define CDUAL_HPP

#include "dual.hpp"

// -----------------------------------------------------------------
// Classes for automatic differentiation of complex functions
// -----------------------------------------------------------------

template <typename Dual>
class CDual {
public:

  Dual real;
  Dual imag;
  // Constructors
  CDual(const Dual &real_, const Dual &imag_)
      : real(real_),
        imag(imag_) {}
};

// addition operators
template <typename Dual>
CDual<Dual> operator+(const CDual<Dual> &dual1, const CDual<Dual> &dual2) {
  return CDual<Dual>(dual1.real + dual2.real, dual1.imag + dual2.imag);
}

template <typename Dual>
CDual<Dual> operator+(const CDual<Dual> &dual, const double &scalar) {
  return CDual<Dual>(dual.real + scalar, dual.imag);
}

template <typename Dual>
CDual<Dual> operator+(const double &scalar, const CDual<Dual> &dual) {
  return dual + scalar;
}

// subtraction operators
template <typename Dual>
CDual<Dual> operator-(const CDual<Dual> &dual1, const CDual<Dual> &dual2) {
  return CDual<Dual>(dual1.real - dual2.real, dual1.imag - dual2.imag);
}

template <typename Dual>
CDual<Dual> operator-(const CDual<Dual> &dual, double scalar) {
  return CDual<Dual>(dual.real - scalar, dual.imag);
}

template <typename Dual>
CDual<Dual> operator-(const double &scalar, const CDual<Dual> &dual) {
  return CDual<Dual>(scalar - dual.real, -1.0 * dual.imag);
}

// multiplication operators
template <typename Dual>
CDual<Dual> operator*(const CDual<Dual> &dual1, const CDual<Dual> &dual2) {
  return CDual<Dual>(dual1.real * dual2.real - dual1.imag * dual2.imag,
                     dual1.real * dual2.imag + dual1.imag * dual2.real);
}

template <typename Dual>
CDual<Dual> operator*(const CDual<Dual> &dual, const double &scalar) {
  return CDual<Dual>(dual.real * scalar, dual.imag * scalar);
}

template <typename Dual>
CDual<Dual> operator*(const double &scalar, const CDual<Dual> &dual) {
  return dual * scalar;
}

// division operators
template <typename Dual>
CDual<Dual> operator/(const CDual<Dual> &dual1, const CDual<Dual> &dual2) {
  const Dual denom = dual2.real * dual2.real + dual2.imag * dual2.imag;
  const Dual renum = dual1.real * dual2.real + dual1.imag * dual2.imag;
  const Dual imnum = dual1.imag * dual2.real - dual1.real * dual2.imag;
  return CDual<Dual>(renum / denom, imnum / denom);
}

template <typename Dual>
CDual<Dual> operator/(const CDual<Dual> &dual, const double &scalar) {
  return CDual<Dual>(dual.real / scalar, dual.imag / scalar);
}

template <typename Dual>
CDual<Dual> operator/(const double &scalar, const CDual<Dual> &dual) {
  const Dual denom = dual.real * dual.real + dual.imag * dual.imag;
  return CDual<Dual>(scalar * dual.real / denom, -scalar * dual.imag / denom);
}

// -----------------------------------------------------------------
// Wrappers for specific derivative calculations
// -----------------------------------------------------------------

using CDual0 = CDual<Dual0>;
using CDual11 = CDual<Dual11>;
using CDual21 = CDual<Dual21>;
using CDual12 = CDual<Dual12>;
using CDual22 = CDual<Dual22>;

#endif
