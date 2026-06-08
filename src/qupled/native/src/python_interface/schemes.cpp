#include "python_interface/schemes.hpp"
#include "esa.hpp"
#include "hf.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "python_interface/util.hpp"
#include "qstls.hpp"
#include "qstlsf0.hpp"
#include "qstlspimc.hpp"
#include "qstlsiet.hpp"
#include "qvsstls.hpp"
#include "qvsstlsf0.hpp"
#include "recstls.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "stlsiet.hpp"
#include "ve.hpp"
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

  int compute() {
    MPIUtil::init();
    const bool status = TScheme::compute();
    isRoot = MPIUtil::isRoot();
    MPIUtil::finalize();
    return status;
  }

  bool isRoot = true;
};

// Type aliases
using PyHF = PyScheme<HF, Input>;
using PyRpa = PyScheme<Rpa, Input>;
using PyESA = PyScheme<ESA, Input>;
using PyRecStls = PyScheme<RecStls, RecStlsInput>;
using PyStls = PyScheme<Stls, StlsInput>;
using PyQstls = PyScheme<Qstls, QstlsInput>;
using PyQstlsF0 = PyScheme<QstlsF0, QstlsPimcInput>;
using PyQstlsPimc = PyScheme<QstlsPimc, QstlsPimcInput>;
using PyStlsIet = PyScheme<StlsIet, StlsIetInput>;
using PyQstlsIet = PyScheme<QstlsIet, QstlsIetInput>;
using PyVSStls = PyScheme<VSStls, VSStlsInput>;
using PyQVSStls = PyScheme<QVSStls, QVSStlsInput>;
using PyQVSStlsF0 = PyScheme<QVSStlsF0, QVSStlsF0Input>;
using PyVE = PyScheme<VE, VEInput>;

// -----------------------------------------------------------------
// Python exposure helpers
// -----------------------------------------------------------------

template <typename T>
py::array getIdr(const T &scheme) {
  return toNdArray2D(scheme.getIdr());
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
py::array getSsfInput(const T &scheme) {
  return toNdArray(scheme.getSsfInput());
}

template <typename T>
py::array getF0Grid(const T &scheme) {
  return toNdArray(scheme.getF0Grid());
}

template <typename T>
py::array getF0Values(const T &scheme) {
  return toNdArray(scheme.getF0Values());
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
py::array getVCoeff(const T &scheme) {
  return toNdArray(scheme.getVCoeff());
}

template <typename T>
py::array getA1Coeff(const T &scheme) {
  return toNdArray(scheme.getA1Coeff());
}

template <typename T>
py::array getMxcL(const T &scheme) {
  return toNdArray(scheme.getMxcL());
}

template <typename T>
py::array getMatsubaraGrid(const T &scheme) {
  return toNdArray(scheme.getMatsubaraGrid());
}

template <typename T>
void exposeBaseSchemeProperties(py::class_<T> &cls) {
  cls.def("compute", &T::compute)
      .def_property_readonly("idr", &getIdr<T>)
      .def_property_readonly("sdr", &getSdr<T>)
      .def_property_readonly("lfc", &getLfc<T>)
      .def_property_readonly("ssf", &getSsf<T>)
      .def_property_readonly("uint", &T::getUInt)
      .def_property_readonly("wvg", &getWvg<T>)
      .def_readonly("is_root", &T::isRoot);
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
void exposeQstlsPimcSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
  cls.def_property_readonly("f0_grid", &getF0Grid<TScheme>)
      .def_property_readonly("f0_values", &getF0Values<TScheme>);
}

template <typename TScheme, typename TInput>
void exposeIetSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
  cls.def_property_readonly("bf", &getBf<TScheme>);
}

template <typename TScheme, typename TInput>
void exposeRecStlsSchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeBaseSchemeProperties(cls);
  cls.def_property_readonly("ssf_input", &getSsfInput<TScheme>);
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

template <typename TScheme, typename TInput>
void exposeVESchemeClass(py::module_ &m, const std::string &className) {
  auto cls =
      py::class_<TScheme>(m, className.c_str()).def(py::init<const TInput &>());
  exposeIterativeSchemeProperties(cls);
  cls.def_property_readonly("free_energy_integrand",
                            &getFreeEnergyIntegrand<TScheme>)
      .def_property_readonly("free_energy_grid", &getFreeEnergyGrid<TScheme>)
      .def_property_readonly("a_coeff", &getVCoeff<TScheme>)
      .def_property_readonly("a1_coeff", &getA1Coeff<TScheme>)
      .def_property_readonly("mxc_l", &getMxcL<TScheme>)
      .def_property_readonly("matsubara_grid", &getMatsubaraGrid<TScheme>)
      .def_property_readonly("a0_coeff", &TScheme::getA0Coeff)
      .def_property_readonly("kxc0", &TScheme::getKxc0)
      .def_property_readonly("mxc_inf", &TScheme::getMxcInf)
      .def_property_readonly("cxc", &TScheme::getCxc)
      .def_property_readonly("omega_m2", &TScheme::getOmegaM2);
}

// -----------------------------------------------------------------
// Entry point
// -----------------------------------------------------------------

namespace pythonWrappers {

  void exposeSchemes(py::module_ &m) {
    exposeBaseSchemeClass<PyHF, Input>(m, "HF");
    exposeBaseSchemeClass<PyRpa, Input>(m, "Rpa");
    exposeBaseSchemeClass<PyESA, Input>(m, "ESA");
    exposeRecStlsSchemeClass<PyRecStls, RecStlsInput>(m, "RecStls");
    exposeIterativeSchemeClass<PyStls, StlsInput>(m, "Stls");
    exposeIterativeSchemeClass<PyQstls, QstlsInput>(m, "Qstls");
    exposeVESchemeClass<PyVE, VEInput>(m, "VE");
    exposeQstlsPimcSchemeClass<PyQstlsF0, QstlsPimcInput>(m, "QstlsF0");
    exposeQstlsPimcSchemeClass<PyQstlsPimc, QstlsPimcInput>(m, "QstlsPimc");
    exposeIetSchemeClass<PyStlsIet, StlsIetInput>(m, "StlsIet");
    exposeIetSchemeClass<PyQstlsIet, QstlsIetInput>(m, "QstlsIet");
    exposeVSSchemeClass<PyVSStls, VSStlsInput>(m, "VSStls");
    exposeVSSchemeClass<PyQVSStls, QVSStlsInput>(m, "QVSStls");
    exposeVSSchemeClass<PyQVSStlsF0, QVSStlsF0Input>(m, "QVSStlsF0");
  }

} // namespace pythonWrappers
