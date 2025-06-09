#ifndef NEW_PYTHON_INTERFACE_UTILITIES_HPP
#define NEW_PYTHON_INTERFACE_UTILITIES_HPP

#include <pybind11/pybind11.h>

namespace NewPythonWrappers {

  // Function to expose utility functions to Python
  void exposeUtilities(pybind11::module_ &m);

} // namespace NewPythonWrappers

#endif
