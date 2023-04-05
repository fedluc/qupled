#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "input.hpp"
#include "stls.hpp"
#include "qstls.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(qupled)
{

  // Wrapper for vector<double>
  class_<std::vector<double>>("vector<double>")
    .def(vector_indexing_suite<std::vector<double>>() );
  
  // Classes to manage the input
  class_<Input>("Input", init<const double, const double, const string>())
    .add_property("coupling",
		  &Input::getCoupling,
		  &Input::setCoupling)
    .add_property("degeneracy",
		  &Input::getDegeneracy,
		  &Input::setDegeneracy)
    .add_property("int2DScheme",
		  &Input::getInt2DScheme,
		  &Input::setInt2DScheme)
    .add_property("threads",
		  &Input::getNThreads,
		  &Input::setNThreads)
    .add_property("theory",
		  &Input::getTheory,
		  &Input::setTheory)
    .def("print", &Input::print)
    .def("isEqual", &Input::isEqual);

  class_<StlsInput, bases<Input>>("StlsInput", init<const double, const double, const string>())
    .add_property("chemicalPotential",
		  &StlsInput::getChemicalPotentialGuess,
		  &StlsInput::setChemicalPotentialGuess)
    .add_property("error",
		  &StlsInput::getErrMin,
		  &StlsInput::setErrMin)
    .add_property("mixing",
		  &StlsInput::getMixingParameter,
		  &StlsInput::setMixingParameter)
    .add_property("iet",
		  &StlsInput::getIETMapping,
		  &StlsInput::setIETMapping)
    .add_property("matsubara",
		  &StlsInput::getNMatsubara,
		  &StlsInput::setNMatsubara)
    .add_property("iterations",
		  &StlsInput::getNIter,
		  &StlsInput::setNIter)
    .add_property("outputFrequency",
		  &StlsInput::getOutIter,
		  &StlsInput::setOutIter)
    .add_property("restartFile",
		  &StlsInput::getRestartFileName,
		  &StlsInput::setRestartFileName)
    .add_property("resolution",
		  &StlsInput::getWaveVectorGridRes,
		  &StlsInput::setWaveVectorGridRes)
    .add_property("cutoff",
		  &StlsInput::getWaveVectorGridCutoff,
		  &StlsInput::setWaveVectorGridCutoff)
    .def("print", &StlsInput::print)
    .def("isEqual", &StlsInput::isEqual);
  
  class_<QstlsInput>("QstlsInput")
    .add_property("setFixedFileName",
		  &QstlsInput::getFixedFileName,
		  &QstlsInput::setFixedFileName)
    .def("print", &QstlsInput::print)
    .def("isEqual", &QstlsInput::isEqual);
    
  // Class to solve classical schemes
  class_<Stls>("Stls", init<const StlsInput>())
    .def("compute", &Stls::compute)
    .add_property("bf", &Stls::getBf)
    .add_property("idr", &Stls::getIdr)
    .add_property("rdf", &Stls::getRdf)
    .add_property("sdr", &Stls::getSdr)
    .add_property("slfc", &Stls::getSlfc)
    .add_property("ssf", &Stls::getSsf)
    .add_property("ssfHF", &Stls::getSsfHF)
    .add_property("uInt", &Stls::getUInt)
    .add_property("wvg", &Stls::getWvg);
    
  // Class to solve quantum schemes
  class_<Qstls, bases<Stls>>("Qstls", init<const StlsInput, const QstlsInput>())
    .def("compute", &Qstls::compute);
  
}
