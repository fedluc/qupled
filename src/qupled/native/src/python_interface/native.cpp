#include "mpi_util.hpp"
#include "python_interface/inputs.hpp"
#include "python_interface/schemes.hpp"
#include "python_interface/utilities.hpp"

#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// --------------------------------------------------------------------
// Initialization and Finalization
// --------------------------------------------------------------------

void qupledInitialization() {
  // Initialize MPI if necessary
  if (!MPIUtil::isInitialized()) { MPIUtil::init(); }
  // Deactivate default GSL error handler
  gsl_set_error_handler_off();
}

void qupledCleanUp() { MPIUtil::finalize(); }

// --------------------------------------------------------------------
// Pybind11 Module Definition
// --------------------------------------------------------------------

PYBIND11_MODULE(native, m) {
  m.doc() = "qupled native Python bindings via pybind11";

  // Initialization
  qupledInitialization();

  // Register finalization
  m.add_object("_cleanup", py::capsule([]() { qupledCleanUp(); }));

  // Bind submodules and classes
  PythonWrappers::exposeInputs(m);
  PythonWrappers::exposeSchemes(m);
  PythonWrappers::exposeUtilities(m);
}