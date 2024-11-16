#ifndef PYTHON_UTIL_HPP
#define PYTHON_UTIL_HPP

#include "vector2D.hpp"
#include "vector3D.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <vector>

// -----------------------------------------------------------------
// Utility functions to convert between Python and C++ arrays
// -----------------------------------------------------------------

namespace pythonUtil {

  namespace bp = boost::python;
  namespace bn = boost::python::numpy;

  // Check if numpy array is stored in row-major order
  void CheckRowMajor(const bn::ndarray &nda);

  // Convert a numpy array to std::vector<double>
  std::vector<double> toVector(const bn::ndarray &nda);

  // Convert a python list to a std::vector<double>
  std::vector<double> toVector(const bp::list &list);

  // Convert a numpy array to Vector2D
  Vector2D toVector2D(const bn::ndarray &nda);

  // Convery a numpy array to std::vector<std::vector<double>>
  std::vector<std::vector<double>> toDoubleVector(const bn::ndarray &nda);

  // Generic converter from vector type to numpy array
  template <typename T>
  bn::ndarray toNdArrayT(const T &v);

  // Convert std::vector<double> to numpy array
  bn::ndarray toNdArray(const std::vector<double> &v);

  // Convert Vector2D to numpy array
  bn::ndarray toNdArray2D(const Vector2D &v);

  // Convert std::vector<std::vector<double>> to numpy array
  bn::ndarray toNdArray2D(const std::vector<std::vector<double>> &v);

  // Convert Vector3D to numpy array
  bn::ndarray toNdArray3D(const Vector3D &v);

} // namespace pythonUtil

#endif
