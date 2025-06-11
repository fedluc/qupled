#ifndef 
PYTHON_INTERFACE_UTILITIES_HPP
#define 
PYTHON_INTERFACE_UTILITIES_HPP

#include <pybind11/pybind11.h>

namespace pythonWrappers {

  // Function to expose utility functions to Python
  void exposeUtilities(pybind11::module_ &m);

} // namespace pythonWrappers

#endif
