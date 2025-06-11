#ifndef PYTHON_INTERFACE_INPUTS_HPP
#define PYTHON_INTERFACE_INPUTS_HPP

#include <pybind11/pybind11.h>

namespace pythonWrappers {

  void exposeInputs(pybind11::module_ &m);

}

#endif