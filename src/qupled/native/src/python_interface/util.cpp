#include "python_interface/util.hpp"
#include "mpi_util.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>

using namespace std;
namespace py = pybind11;

namespace PythonUtil {

  // ---------------------------------------------------------------------
  // Helper: check C-contiguity
  // ---------------------------------------------------------------------
  void CheckRowMajor(const py::array &arr) {
    if (!(arr.flags() & py::array::c_style)) {
      throw std::runtime_error(
          "The numpy array is not stored in row-major (C-contiguous) order.");
    }
  }

  // ---------------------------------------------------------------------
  // Convert 1D numpy array → std::vector<double>
  // ---------------------------------------------------------------------
  vector<double> toVector(const py::array_t<double> &arr) {
    if (arr.ndim() != 1) {
      throw std::runtime_error("Expected 1D numpy array.");
    }
    return vector<double>(arr.data(), arr.data() + arr.size());
  }

  // ---------------------------------------------------------------------
  // Convert Python list → std::vector<double>
  // ---------------------------------------------------------------------
  vector<double> toVector(const py::list &list) {
    vector<double> v(len(list));
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] = py::cast<double>(list[i]);
    }
    return v;
  }

  // ---------------------------------------------------------------------
  // Convert 2D numpy array → Vector2D (row-major assumed)
  // ---------------------------------------------------------------------
  Vector2D toVector2D(const py::array_t<double> &arr) {
    if (arr.ndim() != 2) {
      throw std::runtime_error("Expected 2D numpy array.");
    }
    CheckRowMajor(arr);
    const ssize_t rows = arr.shape(0);
    const ssize_t cols = arr.shape(1);
    Vector2D result(rows, cols);
    const double *data = arr.data();
    for (ssize_t i = 0; i < rows; ++i) {
      for (ssize_t j = 0; j < cols; ++j) {
        result(i, j) = *(data + i * cols + j);
      }
    }
    return result;
  }

  // ---------------------------------------------------------------------
  // Convert 2D numpy array → vector<vector<double>>
  // ---------------------------------------------------------------------
  vector<vector<double>> toDoubleVector(const py::array_t<double> &arr) {
    if (arr.ndim() != 2) {
      throw std::runtime_error("Expected 2D numpy array.");
    }
    CheckRowMajor(arr);

    const ssize_t rows = arr.shape(0);
    const ssize_t cols = arr.shape(1);

    vector<vector<double>> result(rows, vector<double>(cols));
    const double *data = arr.data();

    for (ssize_t i = 0; i < rows; ++i) {
      for (ssize_t j = 0; j < cols; ++j) {
        result[i][j] = *(data + i * cols + j);
      }
    }

    return result;
  }

  // ---------------------------------------------------------------------
  // Convert std::vector<T> → py::array
  // ---------------------------------------------------------------------
  template <typename T>
  py::array_t<double> toNdArrayT(const T &v) {
    return py::array_t<double>(v.size(), v.data());
  }

  py::array toNdArray(const vector<double> &v) { return toNdArrayT(v); }

  // ---------------------------------------------------------------------
  // Convert Vector2D → py::array
  // ---------------------------------------------------------------------
  py::array toNdArray2D(const Vector2D &v) {
    py::array result = toNdArrayT(v);
    result.resize({v.size(0), v.size(1)});
    return result;
  }

  // Convert vector<vector<double>> → py::array via Vector2D
  py::array toNdArray2D(const vector<vector<double>> &v) {
    return toNdArray2D(Vector2D(v));
  }

  // Convert Vector3D → py::array
  py::array toNdArray3D(const Vector3D &v) {
    py::array result = toNdArrayT(v);
    result.resize({v.size(0), v.size(1), v.size(2)});
    return result;
  }

} // namespace PythonUtil
