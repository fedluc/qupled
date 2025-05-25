#include "mpi_util.hpp"
#include "python_interface/inputs.hpp"
#include "python_interface/schemes.hpp"
#include "python_interface/utilities.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <gsl/gsl_errno.h>

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
  // Exposed classes and methods
  pythonWrappers::exposeInputs();
  pythonWrappers::exposeSchemes();
  pythonWrappers::exposeUtilities();
}
