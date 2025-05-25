#include "python_interface/schemes.hpp"
#include "esa.hpp"
#include "hf.hpp"
#include "input.hpp"
#include "python_interface/util.hpp"
#include "qstls.hpp"
#include "qstlsiet.hpp"
#include "qvsstls.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "stlsiet.hpp"
#include "vsstls.hpp"

using namespace pythonUtil;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Template class for Python wrapper
// -----------------------------------------------------------------

template <typename TScheme, typename TInput>
class PyScheme : public TScheme {
public:

  explicit PyScheme(const TInput &in)
      : TScheme(std::make_shared<TInput>(in)) {}
};

using PyHF = PyScheme<HF, Input>;
using PyRpa = PyScheme<Rpa, Input>;
using PyESA = PyScheme<ESA, Input>;
using PyStls = PyScheme<Stls, StlsInput>;
using PyQstls = PyScheme<Qstls, QstlsInput>;
using PyStlsIet = PyScheme<StlsIet, StlsIetInput>;
using PyQstlsIet = PyScheme<QstlsIet, QstlsIetInput>;
using PyVSStls = PyScheme<VSStls, VSStlsInput>;
using PyQVSStls = PyScheme<QVSStls, QVSStlsInput>;

// -----------------------------------------------------------------
// Template functions to expose scheme properties to Python
// -----------------------------------------------------------------

template <typename T>
bn::ndarray getIdr(const T &scheme) {
  return toNdArray2D(scheme.getIdr());
}

template <typename T>
bn::ndarray getRdf(const T &scheme, const bn::ndarray &r) {
  return toNdArray(scheme.getRdf(toVector(r)));
}

template <typename T>
bn::ndarray getSdr(const T &scheme) {
  return toNdArray(scheme.getSdr());
}

template <typename T>
bn::ndarray getLfc(const T &scheme) {
  return toNdArray2D(scheme.getLfc());
}

template <typename T>
bn::ndarray getSsf(const T &scheme) {
  return toNdArray(scheme.getSsf());
}

template <typename T>
bn::ndarray getWvg(const T &scheme) {
  return toNdArray(scheme.getWvg());
}

template <typename T>
bn::ndarray getBf(const T &scheme) {
  return toNdArray(scheme.getBf());
}

template <typename T>
bn::ndarray getFreeEnergyIntegrand(const T &scheme) {
  return toNdArray2D(scheme.getFreeEnergyIntegrand());
}

template <typename T>
bn::ndarray getFreeEnergyGrid(const T &scheme) {
  return toNdArray(scheme.getFreeEnergyGrid());
}

template <typename T>
void exposeBaseSchemeProperties(bp::class_<T> &cls) {
  cls.def("compute", &T::compute);
  cls.def("rdf", &getRdf<T>);
  cls.add_property("idr", &getIdr<T>);
  cls.add_property("sdr", &getSdr<T>);
  cls.add_property("lfc", &getLfc<T>);
  cls.add_property("ssf", &getSsf<T>);
  cls.add_property("uint", &T::getUInt);
  cls.add_property("wvg", &getWvg<T>);
}

template <typename T>
void exposeIterativeSchemeProperties(bp::class_<T> &cls) {
  exposeBaseSchemeProperties(cls);
  cls.add_property("error", &T::getError);
}

template <typename TScheme, typename TInput>
void exposeBaseSchemeClass(const std::string &className) {
  bp::class_<TScheme> cls(className.c_str(), bp::init<const TInput>());
  exposeBaseSchemeProperties(cls);
}

template <typename TScheme, typename TInput>
void exposeIterativeSchemeClass(const std::string &className) {
  bp::class_<TScheme> cls(className.c_str(), bp::init<const TInput>());
  exposeIterativeSchemeProperties(cls);
}

template <typename TScheme, typename TInput>
void exposeIetSchemeClass(const std::string &className) {
  bp::class_<TScheme> cls(className.c_str(), bp::init<const TInput>());
  exposeIterativeSchemeProperties(cls);
  cls.add_property("bf", &getBf<TScheme>);
}

template <typename TScheme, typename TInput>
void exposeVSSchemeClass(const std::string &className) {
  bp::class_<TScheme> cls(className.c_str(), bp::init<const TInput>());
  exposeIterativeSchemeProperties(cls);
  cls.add_property("alpha", &TScheme::getAlpha);
  cls.add_property("free_energy_integrand", &getFreeEnergyIntegrand<TScheme>);
  cls.add_property("free_energy_grid", &getFreeEnergyGrid<TScheme>);
}

// -----------------------------------------------------------------
// All schemes classes exposed to Python
// -----------------------------------------------------------------

namespace pythonWrappers {
  void exposeSchemes() {
    exposeBaseSchemeClass<PyHF, Input>("HF");
    exposeBaseSchemeClass<PyRpa, Input>("Rpa");
    exposeBaseSchemeClass<PyESA, Input>("ESA");
    exposeIterativeSchemeClass<PyStls, StlsInput>("Stls");
    exposeIterativeSchemeClass<PyQstls, QstlsInput>("Qstls");
    exposeIetSchemeClass<PyStlsIet, StlsIetInput>("StlsIet");
    exposeIetSchemeClass<PyQstlsIet, QstlsIetInput>("QstlsIet");
    exposeVSSchemeClass<PyVSStls, VSStlsInput>("VSStls");
    exposeVSSchemeClass<PyQVSStls, QVSStlsInput>("QVSStls");
  }
} // namespace pythonWrappers