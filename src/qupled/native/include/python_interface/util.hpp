#ifndef PYTHON_INTERFACE_UTIL_HPP
#define PYTHON_INTERFACE_UTIL_HPP

#include "vector2D.hpp"
#include "vector3D.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

namespace pythonUtil {

  namespace py = pybind11;
  using py::array;
  using py::array_t;
  using py::list;

  // Check if a numpy array is stored in row-major order
  void CheckRowMajor(const py::array &arr);

  // Convert Python list or 1D array to std::vector<double>
  std::vector<double> toVector(const py::array_t<double> &arr);
  std::vector<double> toVector(const py::list &list);

  // Convert 2D array to Vector2D or std::vector<std::vector<double>>
  Vector2D toVector2D(const py::array_t<double> &arr);
  std::vector<std::vector<double>>
  toDoubleVector(const py::array_t<double> &arr);

  // Convert native C++ containers to numpy arrays
  py::array toNdArray(const std::vector<double> &v);
  py::array toNdArray2D(const Vector2D &v);
  py::array toNdArray2D(const std::vector<std::vector<double>> &v);
  py::array toNdArray3D(const Vector3D &v);

} // namespace pythonUtil

#endif // PYTHON_INTERFACE_UTIL_HPP