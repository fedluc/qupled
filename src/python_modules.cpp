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
    .def("print", &Input::print);


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
    .def("print", &StlsInput::print);
    
  class_<QstlsInput>("QstlsInput")
    .def("setFixedFileName", &QstlsInput::setFixedFileName)
    .def("print", &QstlsInput::print);
    
  // Class to solve classical schemes
  class_<Stls>("Stls", init<const StlsInput>())
    .def("compute", &Stls::compute);

  // Class to solve quantum schemes
  class_<Qstls>("Qstls", init<const StlsInput, const QstlsInput>())
    .def("compute", &Qstls::compute);
  
}
