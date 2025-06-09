#ifndef NEW_PYTHON_INTERFACE_INPUTS_HPP
#define NEW_PYTHON_INTERFACE_INPUTS_HPP

#include <pybind11/pybind11.h>

namespace NewPythonWrappers {

  void exposeInputs(pybind11::module_ &m);

}

#endif