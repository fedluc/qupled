#include "python_inputs.hpp"
#include "python_schemes.hpp"
#include "python_util.hpp"
#include "python_wrappers.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;

// Initialization code for the qupled module
void qupledInitialization() {
  // Initialize MPI if necessary
  if (!MPIUtil::isInitialized()) { MPIUtil::init(); }
  // Deactivate default GSL error handler
  gsl_set_error_handler_off();
}

// Clean up code to call when the python interpreter exists
void qupledCleanUp() { MPIUtil::finalize(); }

// Classes exposed to Python
BOOST_PYTHON_MODULE(native) {

  // Docstring formatting
  bp::docstring_options docopt;
  docopt.enable_all();
  docopt.disable_cpp_signatures();

  // Numpy library initialization
  bn::initialize();

  // Module initialization
  qupledInitialization();

  // Register cleanup function
  std::atexit(qupledCleanUp);

  // Exposes the schemes
  exposeInputs();
  exposeSchemes();

  // MPI class
  bp::class_<PyMPI>("MPI")
      .def("rank", &PyMPI::rank)
      .staticmethod("rank")
      .def("is_root", &PyMPI::isRoot)
      .staticmethod("is_root")
      .def("barrier", &PyMPI::barrier)
      .staticmethod("barrier")
      .def("timer", &PyMPI::timer)
      .staticmethod("timer");

  // Post-process methods
  bp::def("compute_rdf", &PyThermo::computeRdf);
  bp::def("compute_internal_energy", &PyThermo::computeInternalEnergy);
  bp::def("compute_free_energy", &PyThermo::computeFreeEnergy);
}
