#include "python_interface/utilities.hpp"
#include "dimensions_util.hpp"
#include "mpi_util.hpp"
#include "python_interface/util.hpp"
#include "thermo_util.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pythonUtil;

// -----------------------------------------------------------------
// Thermodynamic Utility Functions
// -----------------------------------------------------------------

py::array computeRdf(const py::array_t<double> &rIn,
                     const py::array_t<double> &wvgIn,
                     const py::array_t<double> &ssfIn,
                     const dimensionsUtil::Dimension &dim) {
  const std::vector<double> r = toVector(rIn);
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> ssf = toVector(ssfIn);
  return toNdArray(thermoUtil::computeRdf(r, wvg, ssf, dim));
}

double computeInternalEnergy(const py::array_t<double> &wvgIn,
                             const py::array_t<double> &ssfIn,
                             const double &coupling,
                             const dimensionsUtil::Dimension &dim) {
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> ssf = toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling, dim);
}

double computeFreeEnergy(const py::array_t<double> &gridIn,
                         const py::array_t<double> &rsuIn,
                         const double &coupling) {
  const std::vector<double> grid = toVector(gridIn);
  const std::vector<double> rsu = toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}

// -----------------------------------------------------------------
// All utilities exposed to Python
// -----------------------------------------------------------------

namespace pythonWrappers {

  void exposePostProcessingMethods(py::module_ &m) {
    m.def("compute_rdf", &computeRdf, "Compute radial distribution function");
    m.def("compute_internal_energy",
          &computeInternalEnergy,
          "Compute internal energy");
    m.def("compute_free_energy", &computeFreeEnergy, "Compute free energy");
  }

  void exposeMPIClass(py::module_ &m) {
    m.attr("uses_mpi") = py::bool_(MPIUtil::isUsed);
  }

  void exposeDimensions(py::module_ &m) {
    m.attr("Dimension") = py::enum_<dimensionsUtil::Dimension>(m, "Dimension")
                              .value("D3", dimensionsUtil::Dimension::D3)
                              .value("D2", dimensionsUtil::Dimension::D2);
  }

  void exposeUtilities(py::module_ &m) {
    exposePostProcessingMethods(m);
    exposeMPIClass(m);
    exposeDimensions(m);
  }

} // namespace pythonWrappers
