#ifndef NEW_PYTHON_INTERFACE_SCHEMES_HPP
#define NEW_PYTHON_INTERFACE_SCHEMES_HPP

#include <pybind11/pybind11.h>

namespace NewPythonWrappers {

  // Function to expose schemes to Python
  void exposeSchemes(pybind11::module_ &m);

} // namespace NewPythonWrappers

#endif