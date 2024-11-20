#include "python_util.hpp"
#include "mpi_util.hpp"

using namespace std;
using namespace MPIUtil;

namespace pythonUtil {

  void CheckRowMajor(const bn::ndarray &nda) {
    const bn::ndarray::bitflag flags = nda.get_flags();
    const bool isRowMajor = flags & bn::ndarray::C_CONTIGUOUS;
    if (!isRowMajor) {
      throwError(
          "The numpy array is not stored in row major order (c-contiguous)");
    }
  }

  vector<double> toVector(const bn::ndarray &nda) {
    if (nda.get_nd() != 1) { throwError("Incorrect numpy array dimensions"); }
    const Py_intptr_t *shape = nda.get_shape();
    const int dim = nda.get_nd();
    // the numpy array is flattened to a one dimensional std::vector
    Py_intptr_t n = 1;
    for (int i = 0; i < dim; ++i) {
      n *= shape[i];
    }
    double *ptr = reinterpret_cast<double *>(nda.get_data());
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) {
      v[i] = *(ptr + i);
    }
    return v;
  }

  vector<double> toVector(const bp::list &list) {
    int n = len(list);
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) {
      v[i] = bp::extract<double>(list[i]);
    }
    return v;
  }

  Vector2D toVector2D(const bn::ndarray &nda) {
    if (nda.get_nd() != 2) { throwError("Incorrect numpy array dimensions"); }
    CheckRowMajor(nda);
    const Py_intptr_t *shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    Vector2D v(sz1, sz2);
    double *ptr = reinterpret_cast<double *>(nda.get_data());
    for (int i = 0; i < sz1; ++i) {
      for (int j = 0; j < sz2; ++j) {
        v(i, j) = *(ptr + j + i * sz2);
      }
    }
    return v;
  }

  vector<vector<double>> toDoubleVector(const bn::ndarray &nda) {
    if (nda.get_nd() != 2) { throwError("Incorrect numpy array dimensions"); }
    CheckRowMajor(nda);
    const Py_intptr_t *shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    vector<vector<double>> v(sz1);
    double *ptr = reinterpret_cast<double *>(nda.get_data());
    for (int i = 0; i < sz1; ++i) {
      v[i].resize(sz2);
      for (int j = 0; j < sz2; ++j) {
        v[i][j] = *(ptr + j + i * sz2);
      }
    }
    return v;
  }

  template <typename T>
  bn::ndarray toNdArrayT(const T &v) {
    Py_intptr_t shape[1];
    shape[0] = v.size();
    bn::ndarray result = bn::zeros(1, shape, bn::dtype::get_builtin<double>());
    std::copy(
        v.begin(), v.end(), reinterpret_cast<double *>(result.get_data()));
    return result;
  }

  bn::ndarray toNdArray(const vector<double> &v) { return toNdArrayT(v); }

  bn::ndarray toNdArray2D(const Vector2D &v) {
    bn::ndarray result = toNdArrayT(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1)));
    return result;
  }

  bn::ndarray toNdArray2D(const vector<vector<double>> &v) {
    return toNdArray2D(Vector2D(v));
  }

  bn::ndarray toNdArray3D(const Vector3D &v) {
    bn::ndarray result = toNdArrayT(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1), v.size(2)));
    return result;
  }

} // namespace pythonUtil
