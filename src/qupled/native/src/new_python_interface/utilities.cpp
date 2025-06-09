#include "new_python_interface/utilities.hpp"
#include "mpi_util.hpp"
#include "new_python_interface/util.hpp"
#include "thermo_util.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace NewPythonUtil;

// -----------------------------------------------------------------
// Thermodynamic Utility Functions
// -----------------------------------------------------------------

py::array computeRdf(const py::array_t<double> &rIn,
                     const py::array_t<double> &wvgIn,
                     const py::array_t<double> &ssfIn) {
  const std::vector<double> r = toVector(rIn);
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> ssf = toVector(ssfIn);
  return toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
}

double computeInternalEnergy(const py::array_t<double> &wvgIn,
                             const py::array_t<double> &ssfIn,
                             const double &coupling) {
  const std::vector<double> wvg = toVector(wvgIn);
  const std::vector<double> ssf = toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
}

double computeFreeEnergy(const py::array_t<double> &gridIn,
                         const py::array_t<double> &rsuIn,
                         const double &coupling) {
  const std::vector<double> grid = toVector(gridIn);
  const std::vector<double> rsu = toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}

// -----------------------------------------------------------------
// Static MPI Info Class
// -----------------------------------------------------------------

class PyMPI {
public:

  static int rank() { return MPIUtil::rank(); }
  static bool is_root() { return MPIUtil::isRoot(); }
  static void barrier() { return MPIUtil::barrier(); }
  static double timer() { return MPIUtil::timer(); }
};

// -----------------------------------------------------------------
// All utilities exposed to Python
// -----------------------------------------------------------------

namespace NewPythonWrappers {

  void exposePostProcessingMethods(py::module_ &m) {
    m.def("compute_rdf", &computeRdf, "Compute radial distribution function");
    m.def("compute_internal_energy",
          &computeInternalEnergy,
          "Compute internal energy");
    m.def("compute_free_energy", &computeFreeEnergy, "Compute free energy");
  }

  void exposeMPIClass(py::module_ &m) {
    py::class_<PyMPI>(m, "MPI")
        .def_static("rank", &PyMPI::rank)
        .def_static("is_root", &PyMPI::is_root)
        .def_static("barrier", &PyMPI::barrier)
        .def_static("timer", &PyMPI::timer);
  }

  void exposeUtilities(py::module_ &m) {
    exposePostProcessingMethods(m);
    exposeMPIClass(m);
  }

} // namespace NewPythonWrappers
