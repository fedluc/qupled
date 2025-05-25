#include "python_interface/utilities.hpp"
#include "mpi_util.hpp"
#include "python_interface/util.hpp"
#include "thermo_util.hpp"

using namespace std;
using namespace pythonUtil;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// PyThermo
// -----------------------------------------------------------------

bn::ndarray computeRdf(const bn::ndarray &rIn,
                       const bn::ndarray &wvgIn,
                       const bn::ndarray &ssfIn) {
  const vector<double> &r = toVector(rIn);
  const vector<double> &wvg = toVector(wvgIn);
  const vector<double> &ssf = toVector(ssfIn);
  return toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
}

double computeInternalEnergy(const bn::ndarray &wvgIn,
                             const bn::ndarray &ssfIn,
                             const double &coupling) {
  const vector<double> &wvg = toVector(wvgIn);
  const vector<double> &ssf = toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
}

double computeFreeEnergy(const bn::ndarray &gridIn,
                         const bn::ndarray &rsuIn,
                         const double &coupling) {
  const vector<double> &grid = toVector(gridIn);
  const vector<double> &rsu = toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}

class PyMPI {
public:

  static int rank() { return MPIUtil::rank(); }
  static bool isRoot() { return MPIUtil::isRoot(); }
  static void barrier() { return MPIUtil::barrier(); }
  static double timer() { return MPIUtil::timer(); }
};

// -----------------------------------------------------------------
// All utilities exposed to Python
// -----------------------------------------------------------------

void exposePostProcessingMethods() {
  bp::def("compute_rdf", &computeRdf);
  bp::def("compute_internal_energy", &computeInternalEnergy);
  bp::def("compute_free_energy", &computeFreeEnergy);
}

void exposeMPIClass() {
  bp::class_<PyMPI> cls("MPI");
  cls.def("rank", &PyMPI::rank);
  cls.staticmethod("rank");
  cls.def("is_root", &PyMPI::isRoot);
  cls.staticmethod("is_root");
  cls.def("barrier", &PyMPI::barrier);
  cls.staticmethod("barrier");
  cls.def("timer", &PyMPI::timer);
  cls.staticmethod("timer");
}

namespace pythonWrappers {
  void exposeUtilities() {
    exposePostProcessingMethods();
    exposeMPIClass();
  }
} // namespace pythonWrappers
