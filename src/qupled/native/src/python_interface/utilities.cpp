#include "python_interface/utilities.hpp"
#include "python_interface/util.hpp"
#include "schemes/input.hpp"
#include "thermo/thermo_util.hpp"
#include "util/database.hpp"
#include "util/dimensions_util.hpp"
#include "util/mpi_util.hpp"

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

py::array computeItcfNonInteracting(const Input &in,
                                    const py::array_t<double> &wvgIn,
                                    const py::array_t<double> &tauValuesIn,
                                    const double mu,
                                    const py::array_t<double> &idrIn) {
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> tauValues = toVector(tauValuesIn);
  const Vector2D idr = toVector2D(idrIn);
  return toNdArray2D(thermoUtil::computeItcfNonInteracting(
      std::make_shared<Input>(in), wvg, tauValues, mu, idr));
}

py::array computeItcf(const Input &in,
                      const py::array_t<double> &wvgIn,
                      const py::array_t<double> &tauValuesIn,
                      const double mu,
                      const py::array_t<double> &idrIn,
                      const py::array_t<double> &lfcIn) {
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> tauValues = toVector(tauValuesIn);
  const Vector2D idr = toVector2D(idrIn);
  const Vector2D lfc = toVector2D(lfcIn);
  return toNdArray2D(
      thermoUtil::computeItcf(std::make_shared<Input>(in), wvg, tauValues, mu, idr, lfc));
}

// -----------------------------------------------------------------
// All utilities exposed to Python
// -----------------------------------------------------------------

namespace pythonWrappers {

  void exposePostProcessingMethods(py::module_ &m) {
    m.def("compute_rdf", &computeRdf, "Compute radial distribution function");
    m.def("compute_internal_energy",
          &computeInternalEnergy,
          "Compute the internal energy");
    m.def("compute_free_energy", &computeFreeEnergy, "Compute the free energy");
    m.def("compute_itcf_non_interacting",
          &computeItcfNonInteracting,
          "Compute the non-interacting (Hartree-Fock) imaginary-time correlation function");
    m.def("compute_itcf",
          &computeItcf,
          "Compute the imaginary-time correlation function" );
  }

  void exposeMPIClass(py::module_ &m) {
    m.attr("uses_mpi") = py::bool_(MPIUtil::isUsed);
  }

  void exposeDimensions(py::module_ &m) {
    m.attr("Dimension") = py::enum_<dimensionsUtil::Dimension>(m, "Dimension")
                              .value("D3", dimensionsUtil::Dimension::D3)
                              .value("D2", dimensionsUtil::Dimension::D2);
  }

  void exposeDatabaseMethods(py::module_ &m) {
    m.def("delete_blob_data_on_disk",
          &databaseUtil::deleteBlobDataOnDisk,
          "Delete blob data on disk for a given run ID");
  }

  void exposeUtilities(py::module_ &m) {
    exposePostProcessingMethods(m);
    exposeMPIClass(m);
    exposeDimensions(m);
    exposeDatabaseMethods(m);
  }

} // namespace pythonWrappers
