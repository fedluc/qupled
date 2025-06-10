#ifndef PYTHON_INTERFACE_SCHEMES_HPP
#define PYTHON_INTERFACE_SCHEMES_HPP

#include <pybind11/pybind11.h>

namespace PythonWrappers {

  // Function to expose schemes to Python
  void exposeSchemes(pybind11::module_ &m);

} // namespace PythonWrappers

#endif