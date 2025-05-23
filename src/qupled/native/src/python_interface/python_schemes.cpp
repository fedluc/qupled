#include "python_schemes.hpp"
#include "esa.hpp"
#include "hf.hpp"
#include "input.hpp"
#include "python_util.hpp"
#include "qstls.hpp"
#include "qstlsiet.hpp"
#include "qvsstls.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "stlsiet.hpp"
#include "vsstls.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

using namespace pythonUtil;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Class schemes exposed to Python
// -----------------------------------------------------------------

template <typename Scheme, typename Input>
class PyScheme : public Scheme {
public:

  explicit PyScheme(const Input &in)
      : Scheme(std::make_shared<Input>(in)) {}
};

// -----------------------------------------------------------------
// WORK IN PROGRESS: Start
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

template <typename Scheme>
void exposeBaseSchemeProperties(bp::class_<Scheme> &cls) {
  cls.def("compute", &Scheme::compute)
      .def("rdf", &getRdf<Scheme>)
      .add_property("idr", &getIdr<Scheme>)
      .add_property("sdr", &getSdr<Scheme>)
      .add_property("lfc", &getLfc<Scheme>)
      .add_property("ssf", &getSsf<Scheme>)
      .add_property("uint", &Scheme::getUInt)
      .add_property("wvg", &getWvg<Scheme>);
}

template <typename Scheme>
void exposeIterativeSchemeProperties(bp::class_<Scheme> &cls) {
  exposeBaseSchemeProperties(cls);
  cls.add_property("error", &Scheme::getError);
}

template <typename Scheme, typename Input>
void exposeBaseSchemeClass(const std::string &className) {
  bp::class_<PyScheme<Scheme, Input>> cls(className.c_str(),
                                          bp::init<const Input>());
  exposeBaseSchemeProperties(cls);
}

template <typename Scheme, typename Input>
void exposeIterativeSchemeClass(const std::string &className) {
  bp::class_<PyScheme<Scheme, Input>> cls(className.c_str(),
                                          bp::init<const Input>());
  exposeIterativeSchemeProperties(cls);
}

template <typename Scheme, typename Input>
void exposeIetSchemeClass(const std::string &className) {
  bp::class_<PyScheme<Scheme, Input>> cls(className.c_str(),
                                          bp::init<const Input>());
  exposeIterativeSchemeProperties(cls);
  cls.add_property("bf", &getBf<Scheme>);
}

template <typename Scheme, typename Input>
void exposeVSSchemeClass(const std::string &className) {
  bp::class_<PyScheme<Scheme, Input>> cls(className.c_str(),
                                          bp::init<const Input>());
  exposeIterativeSchemeProperties(cls);
  cls.add_property("alpha", &Scheme::getAlpha);
  cls.add_property("free_energy_integrand", &getFreeEnergyIntegrand<Scheme>);
  cls.add_property("free_energy_grid", &getFreeEnergyGrid<Scheme>);
}

void exposeSchemes() {
  exposeBaseSchemeClass<HF, Input>("HF");
  exposeBaseSchemeClass<Rpa, Input>("Rpa");
  exposeBaseSchemeClass<ESA, Input>("ESA");
  exposeIterativeSchemeClass<Stls, StlsInput>("Stls");
  exposeIterativeSchemeClass<Qstls, QstlsInput>("Qstls");
  exposeIetSchemeClass<StlsIet, StlsIetInput>("StlsIet");
  exposeIetSchemeClass<QstlsIet, QstlsIetInput>("QstlsIet");
  exposeVSSchemeClass<VSStls, VSStlsInput>("VSStls");
  exposeVSSchemeClass<QVSStls, QVSStlsInput>("QVSStls");
}