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
  std::cerr << rpa.getSsf().size() << std::endl;
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

// namespace RpaInputWrapper {
  
//   bn::ndarray getChemicalPotentialGuess(RpaInput &in){
//     return vp::toNdArray(in.getChemicalPotentialGuess());
//   }
  
//   void setChemicalPotentialGuess(RpaInput &in,
// 				 const bp::list &muGuess){
//     in.setChemicalPotentialGuess(arrayWrapper::toVector(muGuess));
//   }

// }

// namespace StlsInputWrapper {
  
//   struct SlfcGuess {
//     bn::ndarray wvg = arrayWrapper::toNdArray(vector<double>(0));
//     bn::ndarray slfc = arrayWrapper::toNdArray(vector<double>(0));
//   };
    
//   StlsInputWrapper::SlfcGuess getGuess(StlsInput &in){
//     StlsInput::SlfcGuess guess_ = in.getGuess();
//     StlsInputWrapper::SlfcGuess guess;
//     guess.wvg = arrayWrapper::toNdArray(guess_.wvg);
//     guess.slfc = arrayWrapper::toNdArray(guess_.slfc);
//     return guess;
//   }
  
//   void setGuess(StlsInput &in,
// 		const StlsInputWrapper::SlfcGuess &guess){
//     StlsInput::SlfcGuess guess_;
//     guess_.wvg = arrayWrapper::toVector(guess.wvg);
//     guess_.slfc = arrayWrapper::toVector(guess.slfc);
//     in.setGuess(guess_);
//   }

// }


// namespace StlsWrapper {

//   bn::ndarray getBf(const Stls &stls){
//     return arrayWrapper::toNdArray(stls.getBf());
//   }

// }


// namespace VSStlsInputWrapper {

//   bn::ndarray getAlphaGuess(VSStlsInput &in){
//     return arrayWrapper::toNdArray(in.getAlphaGuess());
//   }
  
//   void setAlphaGuess(VSStlsInput &in,
// 		     const bp::list &alphaGuess){
//     in.setAlphaGuess(arrayWrapper::toVector(alphaGuess));
//   }
  
//   struct FreeEnergyIntegrand {
//     bn::ndarray grid = arrayWrapper::toNdArray(vector<double>(0));
//     bn::ndarray integrand = arrayWrapper::toNdArray(vector<double>(0));
//   };
  
//   VSStlsInputWrapper::FreeEnergyIntegrand getFreeEnergyIntegrand(VSStlsInput &in){
//     VSStlsInput::FreeEnergyIntegrand fxcIntegrand_ = in.getFreeEnergyIntegrand();
//     VSStlsInputWrapper::FreeEnergyIntegrand fxcIntegrand;
//     fxcIntegrand.grid = arrayWrapper::toNdArray(fxcIntegrand_.grid);
//     fxcIntegrand.integrand = arrayWrapper::toNdArray2D(fxcIntegrand_.integrand);
//     return fxcIntegrand;
//   }
  
//   void setFreeEnergyIntegrand(VSStlsInput &in,
// 			      const VSStlsInputWrapper::FreeEnergyIntegrand &fxcIntegrand){
//     VSStlsInput::FreeEnergyIntegrand fxcIntegrand_;
//     fxcIntegrand_.grid = arrayWrapper::toVector(fxcIntegrand.grid);
//     fxcIntegrand_.integrand = arrayWrapper::toDoubleVector(fxcIntegrand.integrand);
//     in.setFreeEnergyIntegrand(fxcIntegrand_);
//   }
  
// }

// namespace VSStlsWrapper {
    
//   bn::ndarray getFreeEnergyIntegrand(const VSStls &vsstls){
//     return arrayWrapper::toNdArray2D(vsstls.getFreeEnergyIntegrand());
//   }

//   bn::ndarray getFreeEnergyGrid(const VSStls &vsstls){
//     return arrayWrapper::toNdArray(vsstls.getFreeEnergyGrid());
//   }
  
// }

// namespace QstlsInputWrapper {
  
//   struct QstlsGuess {
//     bn::ndarray wvg = arrayWrapper::toNdArray(vector<double>(0));
//     bn::ndarray ssf = arrayWrapper::toNdArray(vector<double>(0));
//     bn::ndarray adr = arrayWrapper::toNdArray(vector<double>(0));
//     int matsubara = 0;
//   };

//   QstlsInputWrapper::QstlsGuess getGuess(QstlsInput &in){
//     QstlsInput::QstlsGuess guess_ = in.getGuess();
//     QstlsInputWrapper::QstlsGuess guess;
//     const int sz1 = guess_.adr.size(0);
//     const int sz2 = guess_.adr.size(1);
//     guess.wvg = arrayWrapper::toNdArray(guess_.wvg);
//     guess.ssf = arrayWrapper::toNdArray(guess_.ssf);
//     bn::ndarray adrTmp = arrayWrapper::toNdArray(guess_.adr);
//     guess.adr = adrTmp.reshape(bp::make_tuple(sz1, sz2));
//     guess.matsubara = guess_.matsubara;
//     return guess;
//   }
  
//   void setGuess(QstlsInput &in,
// 		const QstlsInputWrapper::QstlsGuess &guess){
//     QstlsInput::QstlsGuess guess_;
//     guess_.wvg = arrayWrapper::toVector(guess.wvg);
//     guess_.ssf = arrayWrapper::toVector(guess.ssf);
//     if (guess.adr.shape(0) > 0) {
//       guess_.adr = arrayWrapper::toVector2D(guess.adr);
//     }
//     guess_.matsubara = guess.matsubara;
//     in.setGuess(guess_);
//   }
  
// }

// namespace QstlsWrapper {
    
//   bn::ndarray getAdr(const Qstls &qstls){
//     return arrayWrapper::toNdArray2D(qstls.getAdr());
//   }
  
//   bn::ndarray getAdrFixed(const Qstls &qstls){
//     return arrayWrapper::toNdArray3D(qstls.getAdrFixed());
//   }

// }

// namespace thermoWrapper {

//   bn::ndarray computeRdf(const bn::ndarray &rIn,
// 			 const bn::ndarray &wvgIn,
// 			 const bn::ndarray &ssfIn) {
//     const vector<double> &r = arrayWrapper::toVector(rIn);
//     const vector<double> &wvg = arrayWrapper::toVector(wvgIn);
//     const vector<double> &ssf = arrayWrapper::toVector(ssfIn);
//     return arrayWrapper::toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
//   }

//   double computeInternalEnergy(const bn::ndarray &wvgIn,
// 			       const bn::ndarray &ssfIn,
// 			       const double &coupling) {
//     const vector<double> &wvg = arrayWrapper::toVector(wvgIn);
//     const vector<double> &ssf = arrayWrapper::toVector(ssfIn);
//     return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
//   }

//   double computeFreeEnergy(const bn::ndarray &gridIn,
// 			   const bn::ndarray &rsuIn,
// 			   const double &coupling) {
//     const vector<double> &grid = arrayWrapper::toVector(gridIn);
//     const vector<double> &rsu = arrayWrapper::toVector(rsuIn);
//     return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
//   }
  
// }
