#include "python_interface/schemes.hpp"
#include "esa.hpp"
#include "hf.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "python_interface/util.hpp"
#include "qstls.hpp"
#include "qstlsiet.hpp"
#include "qvsstls.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "stlsiet.hpp"
#include "vsstls.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pythonUtil;

// -----------------------------------------------------------------
// Template wrapper class for Python construction
// -----------------------------------------------------------------

template <typename TScheme, typename TInput>
class PyScheme : public TScheme {
public:

  explicit PyScheme(const TInput &in)
      : TScheme(std::make_shared<TInput>(in)) {}
};

// Type aliases
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
// Python exposure helpers
// -----------------------------------------------------------------

template <typename T>
py::array getIdr(const T &scheme) {
  return toNdArray2D(scheme.getIdr());
}

template <typename T>
py::array getRdf(const T &scheme, const py::array_t<double> &r) {
  return toNdArray(scheme.getRdf(toVector(r)));
}

template <typename T>
py::array getSdr(const T &scheme) {
  return toNdArray(scheme.getSdr());
}

template <typename T>
py::array getLfc(const T &scheme) {
  return toNdArray2D(scheme.getLfc());
}

template <typename T>
py::array getSsf(const T &scheme) {
  return toNdArray(scheme.getSsf());
}

template <typename T>
py::array getWvg(const T &scheme) {
  return toNdArray(scheme.getWvg());
}

template <typename T>
py::array getBf(const T &scheme) {
  return toNdArray(scheme.getBf());
}

template <typename T>
py::array getFreeEnergyIntegrand(const T &scheme) {
  return toNdArray2D(scheme.getFreeEnergyIntegrand());
}

template <typename T>
py::array getFreeEnergyGrid(const T &scheme) {
  return toNdArray(scheme.getFreeEnergyGrid());
}

template <typename T>
void exposeBaseSchemeProperties(py::class_<T> &cls) {
  cls.def("compute", &T::compute)
      .def("rdf", &getRdf<T>)
      .def_property_readonly("idr", &getIdr<T>)
      .def_property_readonly("sdr", &getSdr<T>)
      .def_property_readonly("lfc", &getLfc<T>)
      .def_property_readonly("ssf", &getSsf<T>)
      .def_property_readonly("uint", &T::getUInt)
      .def_property_readonly("wvg", &getWvg<T>)
      .def_property_readonly("is_root", &T::isRunningOnRoot);
}

template <typename T>
void exposeIterativeSchemeProperties(py::class_<T> &cls) {
  exposeBaseSchemeProperties(cls);
  cls.def_property_readonly("error", &T::getError);
}

template <typename TScheme, typename TInput>
void exposeBaseSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeBaseSchemeProperties(cls);
}

template <typename TScheme, typename TInput>
void exposeIterativeSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
}

template <typename TScheme, typename TInput>
void exposeIetSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
  cls.def_property_readonly("bf", &getBf<TScheme>);
}

template <typename TScheme, typename TInput>
void exposeVSSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
  cls.def_property_readonly("alpha", &TScheme::getAlpha)
      .def_property_readonly("free_energy_integrand",
                             &getFreeEnergyIntegrand<TScheme>)
      .def_property_readonly("free_energy_grid", &getFreeEnergyGrid<TScheme>);
}

// -----------------------------------------------------------------
// Entry point
// -----------------------------------------------------------------

namespace pythonWrappers {

  void exposeSchemes(py::module_ &m) {
    exposeBaseSchemeClass<PyHF, Input>(m, "HF");
    exposeBaseSchemeClass<PyRpa, Input>(m, "Rpa");
    exposeBaseSchemeClass<PyESA, Input>(m, "ESA");
    exposeIterativeSchemeClass<PyStls, StlsInput>(m, "Stls");
    exposeIterativeSchemeClass<PyQstls, QstlsInput>(m, "Qstls");
    exposeIetSchemeClass<PyStlsIet, StlsIetInput>(m, "StlsIet");
    exposeIetSchemeClass<PyQstlsIet, QstlsIetInput>(m, "QstlsIet");
    exposeVSSchemeClass<PyVSStls, VSStlsInput>(m, "VSStls");
    exposeVSSchemeClass<PyQVSStls, QVSStlsInput>(m, "QVSStls");
  }

} // namespace pythonWrappers
