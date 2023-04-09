#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "input.hpp"
#include "stls.hpp"
#include "qstls.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;

// Methods that need wrapping to pass arrays between native and python
namespace arrayWrappers {

  vector<double> toVector(const bn::ndarray &nda){
    const Py_intptr_t* shape = nda.get_shape();
    const int dim = nda.get_nd();
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

  template<typename T>
  bn::ndarray toNdArray(const T &v){
    Py_intptr_t shape[1];
    shape[0] = v.size();
    bn::ndarray result = bn::zeros(1, shape, bn::dtype::get_builtin<double>());
    std::copy(v.begin(), v.end(), reinterpret_cast<double*>(result.get_data()));
    return result;
  }
  
  void setChemicalPotentialGuess(StlsInput &in,
				 const bp::list &muGuess){
    in.setChemicalPotentialGuess(toVector(muGuess));
  }
  
  bn::ndarray getBf(const Stls &stls){
    return toNdArray(stls.getBf());
  }

  bn::ndarray getIdr(const Stls &stls){
    const vecUtil::Vector2D &idrNative = stls.getIdr();
    const size_t nx = idrNative.size(0);
    const size_t nl = idrNative.size(1);
    bn::ndarray idr = toNdArray(idrNative);
    bp::tuple shape = bp::make_tuple(nx, nl);
    idr = idr.reshape(shape);
    return idr;
  }

  bn::ndarray getRdf(const Stls &stls,
		     const bn::ndarray &r){
    return toNdArray(stls.getRdf(toVector(r)));
  }
  
  bn::ndarray getSdr(const Stls &stls){
    return toNdArray(stls.getSdr());
  }
  
  bn::ndarray getSlfc(const Stls &stls){
    return toNdArray(stls.getSlfc());
  }
  
  bn::ndarray getSsf(const Stls &stls){
    return toNdArray(stls.getSsf());
  }

  bn::ndarray getSsfHF(const Stls &stls){
    return toNdArray(stls.getSsfHF());
  }
  
  bn::ndarray getWvg(const Stls &stls){
    return toNdArray(stls.getWvg());
  }

  bn::ndarray getAdr(const Qstls &qstls){
    const vecUtil::Vector2D &adrNative = qstls.getAdr();
    const size_t nx = adrNative.size(0);
    const size_t nl = adrNative.size(1);
    bn::ndarray adr = toNdArray(adrNative);
    bp::tuple shape = bp::make_tuple(nx, nl);
    adr = adr.reshape(shape);
    return adr;
  }
  
  bn::ndarray getAdrFixed(const Qstls &qstls){
    const vecUtil::Vector3D &adrNative = qstls.getAdrFixed();
    const size_t nx = adrNative.size(0);
    const size_t nl = adrNative.size(1);
    bn::ndarray adr = toNdArray(adrNative);
    bp::tuple shape = bp::make_tuple(nx, nl, nx);
    adr = adr.reshape(shape);
    return adr;
  }

}



// Classes exposed to Python
BOOST_PYTHON_MODULE(qupled)
{

  // Numpy library initialization
  bn::initialize();
    
  // Wrapper for vector<double>
  bp::class_<std::vector<double>>("vector<double>")
    .def(bp::vector_indexing_suite<std::vector<double>>() );
  
  // Classes to manage the input
  bp::class_<Input>("Input",
		    bp::init<const double, const double, const string>())
    .add_property("coupling",
		  &Input::getCoupling,
		  &Input::setCoupling)
    .add_property("degeneracy",
		  &Input::getDegeneracy,
		  &Input::setDegeneracy)
    .add_property("int2DScheme",
		  &Input::getInt2DScheme,
		  &Input::setInt2DScheme)
    .add_property("threads",
		  &Input::getNThreads,
		  &Input::setNThreads)
    .add_property("theory",
		  &Input::getTheory,
		  &Input::setTheory)
    .def("print", &Input::print)
    .def("isEqual", &Input::isEqual);

  bp::class_<StlsInput, bp::bases<Input>>("StlsInput",
					  bp::init<const double, const double, const string>())
    .add_property("chemicalPotential",
		  &StlsInput::getChemicalPotentialGuess,
		  &arrayWrappers::setChemicalPotentialGuess)
    .add_property("error",
		  &StlsInput::getErrMin,
		  &StlsInput::setErrMin)
    .add_property("mixing",
		  &StlsInput::getMixingParameter,
		  &StlsInput::setMixingParameter)
    .add_property("iet",
		  &StlsInput::getIETMapping,
		  &StlsInput::setIETMapping)
    .add_property("matsubara",
		  &StlsInput::getNMatsubara,
		  &StlsInput::setNMatsubara)
    .add_property("iterations",
		  &StlsInput::getNIter,
		  &StlsInput::setNIter)
    .add_property("outputFrequency",
		  &StlsInput::getOutIter,
		  &StlsInput::setOutIter)
    .add_property("restartFile",
		  &StlsInput::getRestartFileName,
		  &StlsInput::setRestartFileName)
    .add_property("resolution",
		  &StlsInput::getWaveVectorGridRes,
		  &StlsInput::setWaveVectorGridRes)
    .add_property("cutoff",
		  &StlsInput::getWaveVectorGridCutoff,
		  &StlsInput::setWaveVectorGridCutoff)
    .def("print", &StlsInput::print)
    .def("isEqual", &StlsInput::isEqual);
  
  bp::class_<QstlsInput>("QstlsInput")
    .add_property("setFixedFileName",
		  &QstlsInput::getFixedFileName,
		  &QstlsInput::setFixedFileName)
    .def("print", &QstlsInput::print)
    .def("isEqual", &QstlsInput::isEqual);
    
  // Class to solve classical schemes
  bp::class_<Stls>("Stls",
		   bp::init<const StlsInput>())
    .def("compute", &Stls::compute)
    .def("getRdf", arrayWrappers::getRdf)
    .add_property("bf", arrayWrappers::getBf)
    .add_property("idr", arrayWrappers::getIdr)
    .add_property("sdr", arrayWrappers::getSdr)
    .add_property("slfc", arrayWrappers::getSlfc)
    .add_property("ssf", arrayWrappers::getSsf)
    .add_property("ssfHF", arrayWrappers::getSsfHF)
    .add_property("uInt", &Stls::getUInt)
    .add_property("wvg", arrayWrappers::getWvg);
  
  // Class to solve quantum schemes
  bp::class_<Qstls, bp::bases<Stls>>("Qstls",
				     bp::init<const StlsInput, const QstlsInput>())
    .def("compute", &Qstls::compute)
    .add_property("adr", arrayWrappers::getAdr)
    .add_property("adr_fixed", arrayWrappers::getAdrFixed);
  
}
