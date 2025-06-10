#ifndef PYTHON_INTERFACE_SCHEMES_HPP
#define PYTHON_INTERFACE_SCHEMES_HPP

#include <pybind11/pybind11.h>

namespace pythonWrappers {

  // Function to expose schemes to Python
  void exposeSchemes(pybind11::module_ &m);

} // namespace pythonWrappers

#endif