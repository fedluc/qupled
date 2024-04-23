#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "util.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "python_wrappers.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
namespace vp = vecUtil::python;

// Initialization code for the qupled module
void qupledInitialization() {
  // Initialize MPI if necessary
  if (!parallelUtil::MPI::isInitialized()) {
    parallelUtil::MPI::init();
  }
  // Deactivate default GSL error handler
  gsl_set_error_handler_off();
}

// Clean up code to call when the python interpreter exists
void qupledCleanUp() {
  parallelUtil::MPI::finalize();
}

// Classes exposed to Python
BOOST_PYTHON_MODULE(qupled)
{

  // Docstring formatting
  bp::docstring_options docopt;
  docopt.enable_all();
  docopt.disable_cpp_signatures();
	
  // Numpy library initialization
  bn::initialize();

  // Module initialization
  qupledInitialization();

  // Register cleanup function
  std::atexit(qupledCleanUp);
    
  // Class for the input of the Rpa scheme
  bp::class_<RpaInput>("RpaInput")
    .add_property("coupling",
		  &RpaInput::getCoupling,
		  &RpaInput::setCoupling)
    .add_property("degeneracy",
		  &RpaInput::getDegeneracy,
		  &RpaInput::setDegeneracy)
    .add_property("int2DScheme",
		  &RpaInput::getInt2DScheme,
		  &RpaInput::setInt2DScheme)
    .add_property("intError",
		  &RpaInput::getIntError,
		  &RpaInput::setIntError)
    .add_property("threads",
		  &RpaInput::getNThreads,
		  &RpaInput::setNThreads)
    .add_property("theory",
		  &RpaInput::getTheory,
		  &RpaInput::setTheory)
    .add_property("chemicalPotential",
		  &PyRpaInput::getChemicalPotentialGuess,
		  &PyRpaInput::setChemicalPotentialGuess)
    .add_property("matsubara",
		  &RpaInput::getNMatsubara,
		  &RpaInput::setNMatsubara)
    .add_property("resolution",
		  &RpaInput::getWaveVectorGridRes,
		  &RpaInput::setWaveVectorGridRes)
    .add_property("cutoff",
		  &RpaInput::getWaveVectorGridCutoff,
		  &RpaInput::setWaveVectorGridCutoff)
    .def("print", &RpaInput::print)
    .def("isEqual", &RpaInput::isEqual);

  // Class for the initial guess of the Stls scheme
  bp::class_<StlsInput::SlfcGuess>("SlfcGuess")
    .add_property("wvg",
		  &PySlfcGuess::getWvg,
		  &PySlfcGuess::setWvg)
    .add_property("slfc",
		  &PySlfcGuess::getSlfc,
		  &PySlfcGuess::setSlfc);

  // Class for the input of the Stls scheme
  bp::class_<StlsInput, bp::bases<RpaInput>>("StlsInput")
    .add_property("error",
		  &StlsInput::getErrMin,
		  &StlsInput::setErrMin)
    .add_property("mixing",
		  &StlsInput::getMixingParameter,
		  &StlsInput::setMixingParameter)
    .add_property("iet",
		  &StlsInput::getIETMapping,
		  &StlsInput::setIETMapping)
    .add_property("iterations",
		  &StlsInput::getNIter,
		  &StlsInput::setNIter)
    .add_property("outputFrequency",
		  &StlsInput::getOutIter,
		  &StlsInput::setOutIter)
    .add_property("recoveryFile",
		  &StlsInput::getRecoveryFileName,
		  &StlsInput::setRecoveryFileName)
    .add_property("guess",
		  &StlsInput::getGuess,
		  &StlsInput::setGuess)
    .def("print", &StlsInput::print)
    .def("isEqual", &StlsInput::isEqual);

  // Class for the free energy integrand of the VSStls scheme
  bp::class_<VSStlsInput::FreeEnergyIntegrand>("FreeEnergyIntegrand")
    .add_property("grid",
		  &PyFreeEnergyIntegrand::getGrid,
		  &PyFreeEnergyIntegrand::setGrid)
    .add_property("integrand",
		  &PyFreeEnergyIntegrand::getIntegrand,
		  &PyFreeEnergyIntegrand::setIntegrand);

  // Class for the input of the VSStls scheme
  bp::class_<VSStlsInput, bp::bases<StlsInput>>("VSStlsInput")
    .add_property("errorAlpha",
		  &VSStlsInput::getErrMinAlpha,
		  &VSStlsInput::setErrMinAlpha)
    .add_property("iterationsAlpha",
		  &VSStlsInput::getNIterAlpha,
		  &VSStlsInput::setNIterAlpha)
    .add_property("alpha",
		  &PyVSStlsInput::getAlphaGuess,
		  &PyVSStlsInput::setAlphaGuess)
    .add_property("couplingResolution",
		  &VSStlsInput::getCouplingResolution,
		  &VSStlsInput::setCouplingResolution)
    .add_property("degeneracyResolution",
		  &VSStlsInput::getDegeneracyResolution,
		  &VSStlsInput::setDegeneracyResolution)
    .add_property("freeEnergyIntegrand",
		  &VSStlsInput::getFreeEnergyIntegrand,
		  &VSStlsInput::setFreeEnergyIntegrand)
    .def("print", &VSStlsInput::print)
    .def("isEqual", &VSStlsInput::isEqual);

  // Class for the initial guess of the Qstls scheme
  bp::class_<QstlsInput::QstlsGuess>("QstlsGuess")
    .add_property("wvg",
		  &PyQstlsGuess::getWvg,
		  &PyQstlsGuess::setWvg)
    .add_property("ssf",
		  &PyQstlsGuess::getSsf,
		  &PyQstlsGuess::setSsf)
    .add_property("adr",
		  &PyQstlsGuess::getAdr,
		  &PyQstlsGuess::setAdr)
    .add_property("matsubara",
		  &PyQstlsGuess::getMatsubara,
		  &PyQstlsGuess::setMatsubara);
    
  // Class for the input of the Qstls scheme
  bp::class_<QstlsInput, bp::bases<StlsInput>>("QstlsInput")
    .add_property("guess",
		  &QstlsInput::getGuess,
		  &QstlsInput::setGuess)
    .add_property("fixed",
		  &QstlsInput::getFixed,
		  &QstlsInput::setFixed)
    .add_property("fixediet",
		  &QstlsInput::getFixedIet,
		  &QstlsInput::setFixedIet)
    .def("print", &QstlsInput::print)
    .def("isEqual", &QstlsInput::isEqual);

// Class for the input of the QVSStls scheme
  bp::class_<QVSStlsInput, bp::bases<QstlsInput>>("QVSStlsInput")
    .add_property("errorAlpha",
		  &QVSStlsInput::getErrMinAlpha,
		  &QVSStlsInput::setErrMinAlpha)
    .add_property("iterationsAlpha",
		  &QVSStlsInput::getNIterAlpha,
		  &QVSStlsInput::setNIterAlpha)
    .add_property("alpha",
		  &PyQVSStlsInput::getAlphaGuess,
		  &PyQVSStlsInput::setAlphaGuess)
    .add_property("couplingResolution",
		  &QVSStlsInput::getCouplingResolution,
		  &QVSStlsInput::setCouplingResolution)
    .add_property("degeneracyResolution",
		  &QVSStlsInput::getDegeneracyResolution,
		  &QVSStlsInput::setDegeneracyResolution)
    .add_property("freeEnergyIntegrand",
		  &QVSStlsInput::getFreeEnergyIntegrand,
		  &QVSStlsInput::setFreeEnergyIntegrand)
    .def("print", &QVSStlsInput::print)
    .def("isEqual", &QVSStlsInput::isEqual);

  // Class to solve the classical RPA scheme
  bp::class_<Rpa>("Rpa",
		  bp::init<const RpaInput>())
    .def("compute", &PyRpa::compute)
    .def("rdf", &PyRpa::getRdf)
    .add_property("idr", &PyRpa::getIdr)
    .add_property("sdr", &PyRpa::getSdr)
    .add_property("slfc", &PyRpa::getSlfc)
    .add_property("ssf", &PyRpa::getSsf)
    .add_property("ssfHF", &PyRpa::getSsfHF)
    .add_property("uInt", &PyRpa::getUInt)
    .add_property("wvg", &PyRpa::getWvg)
    .add_property("recovery", &PyRpa::getRecoveryFileName);

  // Class to solve the classical ESA scheme
  bp::class_<ESA, bp::bases<Rpa>>("ESA",
				  bp::init<const RpaInput>())
    .def("compute", &ESA::compute);

  // Class to solve classical schemes
  bp::class_<Stls, bp::bases<Rpa>>("Stls",
				   bp::init<const StlsInput>())
    .def("compute", &PyStls::compute)
    .add_property("error", &PyStls::getError)
    .add_property("bf", &PyStls::getBf);

  // Class to solve the classical VS scheme
  bp::class_<VSStls, bp::bases<Rpa>>("VSStls",
				     bp::init<const VSStlsInput>())
    .def("compute", &PyVSStls::compute)
    .add_property("error", &PyVSStls::getError)
    .add_property("Alpha", &PyVSStls::getAlpha)
    .add_property("freeEnergyIntegrand", &PyVSStls::getFreeEnergyIntegrand)
    .add_property("freeEnergyGrid", &PyVSStls::getFreeEnergyGrid);
      
  // Class to solve quantum schemes
  bp::class_<Qstls, bp::bases<Stls>>("Qstls",
				     bp::init<const QstlsInput>())
    .def("compute", &PyQstls::compute)
    .add_property("error", &PyQstls::getError)
    .add_property("adr", &PyQstls::getAdr);

  // Class to solve the quantum VS scheme
  bp::class_<QVSStls, bp::bases<Rpa>>("QVSStls",
				      bp::init<const QVSStlsInput>())
    .def("compute", &PyQVSStls::compute)
    .add_property("error", &PyQVSStls::getError)
    .add_property("freeEnergyIntegrand", &PyQVSStls::getFreeEnergyIntegrand)
    .add_property("freeEnergyGrid", &PyQVSStls::getFreeEnergyGrid)
    .add_property("adr", &PyQVSStls::getAdr)
    .add_property("Alpha", &PyQVSStls::getAlpha);

  // MPI class
  bp::class_<PyMPI>("MPI")
    .def("rank", &PyMPI::rank)
    .def("isRoot", &PyMPI::isRoot)
    .def("barrier", &PyMPI::barrier)
    .def("timer", &PyMPI::timer);
  
  
  // Post-process methods
  bp::def("computeRdf", &PyThermo::computeRdf);
  bp::def("computeInternalEnergy", &PyThermo::computeInternalEnergy);
  bp::def("computeFreeEnergy", &PyThermo::computeFreeEnergy);

  
}

