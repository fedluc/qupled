#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "util.hpp"
#include "python_wrappers.hpp"

using namespace thermoUtil;
// Methods that need wrapping to pass arrays between native and python

namespace vp {

  void CheckRowMajor(const bn::ndarray &nda) {
    const bn::ndarray::bitflag flags = nda.get_flags();
    const bool isRowMajor = flags & bn::ndarray::C_CONTIGUOUS;
    if (!isRowMajor) {
      throw runtime_error("The numpy array is not stored in row major order (c-contiguous)");
    }
  }
  
  vector<double> toVector(const bn::ndarray &nda){
    if (nda.get_nd() != 1) {
      throw runtime_error("Incorrect numpy array dimensions");
    }
    const Py_intptr_t* shape = nda.get_shape();
    const int dim = nda.get_nd();
    // the numpy array is flattened to a one dimensional std::vector
    Py_intptr_t n = 1;
    for (int i = 0; i < dim; ++i){ n *= shape[i]; }
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) { v[i] = *(ptr + i); }
    return v;
  }

  vector<double> toVector(const bp::list &list){
    int n = len(list);
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i){ v[i] = bp::extract<double>(list[i]); }
    return v;
  }

  vecUtil::Vector2D toVector2D(const bn::ndarray &nda){
    if (nda.get_nd() != 2) {
      throw runtime_error("Incorrect numpy array dimensions");
    }
    CheckRowMajor(nda);
    const Py_intptr_t* shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    vecUtil::Vector2D v(sz1, sz2);
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    for (int i = 0; i < sz1; ++i){
      for (int j = 0; j < sz2; ++j) {
	v(i,j) = *(ptr + j + i*sz2);
      }
    }
    return v;
  }

  vector<vector<double>> toDoubleVector(const bn::ndarray &nda){
    if (nda.get_nd() != 2) {
      throw runtime_error("Incorrect numpy array dimensions");
    }
    CheckRowMajor(nda);
    const Py_intptr_t* shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    vector<vector<double>> v(sz1);
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    for (int i = 0; i < sz1; ++i){
      v[i].resize(sz2);
      for (int j = 0; j < sz2; ++j) {
	v[i][j] = *(ptr + j + i*sz2);
      }
    }
    return v;
  }
  
  template<typename T>
  bn::ndarray toNdArray(const T &v){
    Py_intptr_t shape[1];
    shape[0] = v.size();
    bn::ndarray result = bn::zeros(1, shape, bn::dtype::get_builtin<double>());
    std::copy(v.begin(), v.end(), reinterpret_cast<double*>(result.get_data()));
    return result;
  }
  
  bn::ndarray toNdArray2D(const vecUtil::Vector2D &v){
    bn::ndarray result = toNdArray(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1)));
    return result;
  }

  bn::ndarray toNdArray2D(const vector<vector<double>> &v){
    return toNdArray2D(vecUtil::Vector2D(v));
  }
  
  bn::ndarray toNdArray3D(const vecUtil::Vector3D &v){
    bn::ndarray result = toNdArray(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1), v.size(2)));
    return result;
  }
  
}

// -----------------------------------------------------------------
// PyRpa
// -----------------------------------------------------------------

bn::ndarray PyRpa::getIdr(const Rpa& rpa) {
  return vp::toNdArray2D(rpa.getIdr());
}

bn::ndarray PyRpa::getRdf(const Rpa& rpa,
			  const bn::ndarray &r) {
  return vp::toNdArray(rpa.getRdf(vp::toVector(r)));
}

bn::ndarray PyRpa::getSdr(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSdr());
}

bn::ndarray PyRpa::getSlfc(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSlfc());
}

bn::ndarray PyRpa::getSsf(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSsf());
}

bn::ndarray PyRpa::getSsfHF(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSsfHF());
}

bn::ndarray PyRpa::getWvg(const Rpa& rpa) {
  return vp::toNdArray(rpa.getWvg());
}

double PyRpa::getUInt(const Rpa& rpa) {
  return rpa.getUInt();
}

string PyRpa::getRecoveryFileName(const Rpa& rpa) {
  return rpa.getRecoveryFileName();
}

// -----------------------------------------------------------------
// PyStls
// -----------------------------------------------------------------

int PyStls::compute(Stls& stls) {
  return stls.compute();
}

bn::ndarray PyStls::getBf(const Stls& stls) {
  return vp::toNdArray(stls.getBf());
}

// -----------------------------------------------------------------
// PyVSStls
// -----------------------------------------------------------------

int PyVSStls::compute(VSStls& vsstls) {
  return vsstls.compute();
}

bn::ndarray PyVSStls::getFreeEnergyIntegrand(const VSStls &vsstls){
  return vp::toNdArray2D(vsstls.getFreeEnergyIntegrand());
}

bn::ndarray PyVSStls::getFreeEnergyGrid(const VSStls &vsstls){
  return vp::toNdArray(vsstls.getFreeEnergyGrid());
}

// -----------------------------------------------------------------
// PyQstls
// -----------------------------------------------------------------

int PyQstls::compute(Qstls& qstls) {
  return qstls.compute();
}

bn::ndarray PyQstls::getAdr(const Qstls& qstls) {
  return vp::toNdArray2D(qstls.getAdr());
}

// -----------------------------------------------------------------
// PyThermo
// -----------------------------------------------------------------

bn::ndarray PyThermo::computeRdf(const bn::ndarray &rIn,
				 const bn::ndarray &wvgIn,
				 const bn::ndarray &ssfIn) {
  const vector<double> &r = vp::toVector(rIn);
  const vector<double> &wvg = vp::toVector(wvgIn);
  const vector<double> &ssf = vp::toVector(ssfIn);
  return vp::toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
}

double PyThermo::computeInternalEnergy(const bn::ndarray &wvgIn,
				       const bn::ndarray &ssfIn,
				       const double &coupling) {
  const vector<double> &wvg = vp::toVector(wvgIn);
  const vector<double> &ssf = vp::toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
}

double PyThermo::computeFreeEnergy(const bn::ndarray &gridIn,
				   const bn::ndarray &rsuIn,
				   const double &coupling) {
  const vector<double> &grid = vp::toVector(gridIn);
  const vector<double> &rsu = vp::toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}

