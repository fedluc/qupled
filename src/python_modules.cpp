#include <boost/python.hpp>
#include "input.hpp"
#include "stls.hpp"
#include "qstls.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(qupled)
{

  // Classes to manage the input
  class_<Input>("Input", init<const double, const double, const string>())
    .def("setCoupling", &Input::setCoupling)
    .def("setDegeneracy", &Input::setDegeneracy)
    .def("setInt2DScheme", &Input::setInt2DScheme)
    .def("setNThreads", &Input::setNThreads)
    .def("setTheory", &Input::setTheory)
    .def("getCoupling", &Input::getCoupling)
    .def("getDegeneracy", &Input::getDegeneracy)
    .def("getInt2DScheme", &Input::getInt2DScheme)
    .def("getNThreads", &Input::getNThreads)
    .def("getTheory", &Input::getTheory)
    .def("print", &Input::print)
    .def("isEqual", &Input::isEqual);

  class_<StlsInput, bases<Input>>("StlsInput", init<const double, const double, const string>())
    .def("setChemicalPotentialGuess", &StlsInput::setChemicalPotentialGuess)
    .def("setErrMin", &StlsInput::setErrMin)
    .def("setMixingParameter", &StlsInput::setMixingParameter)
    .def("setIETMapping", &StlsInput::setIETMapping)
    .def("setNMatsubara", &StlsInput::setNMatsubara)
    .def("setNIter", &StlsInput::setNIter)
    .def("setOutIter", &StlsInput::setOutIter)
    .def("setRestartFileName", &StlsInput::setRestartFileName)
    .def("setWaveVectorGridRes", &StlsInput::setWaveVectorGridRes)
    .def("setWaveVectorGridCutoff", &StlsInput::setWaveVectorGridCutoff)
    .def("getChemicalPotentialGuess", &StlsInput::getChemicalPotentialGuess)
    .def("getErrMin", &StlsInput::getErrMin)
    .def("getMixingParameter", &StlsInput::getMixingParameter)
    .def("getIETMapping", &StlsInput::getIETMapping)
    .def("getNMatsubara", &StlsInput::getNMatsubara)
    .def("getNIter", &StlsInput::getNIter)
    .def("getOutIter", &StlsInput::getOutIter)
    .def("getRestartFileName", &StlsInput::getRestartFileName)
    .def("getWaveVectorGridRes", &StlsInput::getWaveVectorGridRes)
    .def("getWaveVectorGridCutoff", &StlsInput::getWaveVectorGridCutoff)
    .def("print", &StlsInput::print)
    .def("isEqual", &StlsInput::isEqual);
  
  class_<QstlsInput>("QstlsInput")
    .def("setFixedFileName", &QstlsInput::setFixedFileName)
    .def("getFixedFileName", &QstlsInput::getFixedFileName)
    .def("print", &QstlsInput::print)
    .def("isEqual", &QstlsInput::isEqual);
    
  // Class to solve classical schemes
  class_<Stls>("Stls", init<const StlsInput>())
    .def("compute", &Stls::compute);

  // Class to solve quantum schemes
  class_<Qstls>("Qstls", init<const StlsInput, const QstlsInput>())
    .def("compute", &Qstls::compute);
  
}
