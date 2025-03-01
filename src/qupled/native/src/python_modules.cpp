#include "mpi_util.hpp"
#include "python_wrappers.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;

// Initialization code for the qupled module
void qupledInitialization() {
  // Initialize MPI if necessary
  if (!MPIUtil::isInitialized()) { MPIUtil::init(); }
  // Deactivate default GSL error handler
  gsl_set_error_handler_off();
}

// Clean up code to call when the python interpreter exists
void qupledCleanUp() { MPIUtil::finalize(); }

// Classes exposed to Python
BOOST_PYTHON_MODULE(native) {

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
      .add_property("coupling", &RpaInput::getCoupling, &RpaInput::setCoupling)
      .add_property(
          "degeneracy", &RpaInput::getDegeneracy, &RpaInput::setDegeneracy)
      .add_property("integral_strategy",
                    &RpaInput::getInt2DScheme,
                    &RpaInput::setInt2DScheme)
      .add_property(
          "integral_error", &RpaInput::getIntError, &RpaInput::setIntError)
      .add_property("threads", &RpaInput::getNThreads, &RpaInput::setNThreads)
      .add_property("theory", &RpaInput::getTheory, &RpaInput::setTheory)
      .add_property("chemical_potential",
                    &PyRpaInput::getChemicalPotentialGuess,
                    &PyRpaInput::setChemicalPotentialGuess)
      .add_property(
          "matsubara", &RpaInput::getNMatsubara, &RpaInput::setNMatsubara)
      .add_property("resolution",
                    &RpaInput::getWaveVectorGridRes,
                    &RpaInput::setWaveVectorGridRes)
      .add_property("cutoff",
                    &RpaInput::getWaveVectorGridCutoff,
                    &RpaInput::setWaveVectorGridCutoff)
      .add_property("frequency_cutoff",
                    &RpaInput::getFrequencyCutoff,
                    &RpaInput::setFrequencyCutoff)
      .def("print", &RpaInput::print)
      .def("is_equal", &RpaInput::isEqual);

  // Class for the input of the Stls scheme
  bp::class_<StlsInput, bp::bases<RpaInput>>("StlsInput")
      .add_property("error", &StlsInput::getErrMin, &StlsInput::setErrMin)
      .add_property("guess", &StlsInput::getGuess, &StlsInput::setGuess)
      .add_property(
          "mapping", &StlsInput::getIETMapping, &StlsInput::setIETMapping)
      .add_property("mixing",
                    &StlsInput::getMixingParameter,
                    &StlsInput::setMixingParameter)
      .add_property("iterations", &StlsInput::getNIter, &StlsInput::setNIter)
      .add_property(
          "output_frequency", &StlsInput::getOutIter, &StlsInput::setOutIter)
      .add_property("recovery_file",
                    &StlsInput::getRecoveryFileName,
                    &StlsInput::setRecoveryFileName)
      .def("print", &StlsInput::print)
      .def("is_equal", &StlsInput::isEqual);

  // Class for the input of the VSStls scheme
  bp::class_<VSInput>("VSInput")
      .add_property(
          "error_alpha", &VSInput::getErrMinAlpha, &VSInput::setErrMinAlpha)
      .add_property(
          "iterations_alpha", &VSInput::getNIterAlpha, &VSInput::setNIterAlpha)
      .add_property(
          "alpha", &PyVSInput::getAlphaGuess, &PyVSInput::setAlphaGuess)
      .add_property("coupling_resolution",
                    &VSInput::getCouplingResolution,
                    &VSInput::setCouplingResolution)
      .add_property("degeneracy_resolution",
                    &VSInput::getDegeneracyResolution,
                    &VSInput::setDegeneracyResolution)
      .add_property("free_energy_integrand",
                    &VSInput::getFreeEnergyIntegrand,
                    &VSInput::setFreeEnergyIntegrand)
      .def("print", &VSInput::print)
      .def("is_equal", &VSInput::isEqual);

  // Class for the input of the VSStls scheme
  bp::class_<VSStlsInput, bp::bases<VSInput, StlsInput>>("VSStlsInput")
      .def("print", &VSStlsInput::print)
      .def("is_equal", &VSStlsInput::isEqual);

  // Class for the input of the Qstls scheme
  bp::class_<QstlsInput, bp::bases<StlsInput>>("QstlsInput")
      .add_property("guess", &QstlsInput::getGuess, &QstlsInput::setGuess)
      .add_property("fixed", &QstlsInput::getFixed, &QstlsInput::setFixed)
      .add_property(
          "fixed_iet", &QstlsInput::getFixedIet, &QstlsInput::setFixedIet)
      .def("print", &QstlsInput::print)
      .def("is_equal", &QstlsInput::isEqual);

  // Class for the input of the QVSStls scheme
  bp::class_<QVSStlsInput, bp::bases<VSInput, QstlsInput>>("QVSStlsInput")
      .def("print", &QVSStlsInput::print)
      .def("is_equal", &QVSStlsInput::isEqual);

  // Class for the initial guess of the Stls scheme
  bp::class_<StlsInput::Guess>("StlsGuess")
      .add_property("wvg", &PyStlsGuess::getWvg, &PyStlsGuess::setWvg)
      .add_property("slfc", &PyStlsGuess::getSlfc, &PyStlsGuess::setSlfc);

  // Class for the initial guess of the Qstls scheme
  bp::class_<QstlsInput::Guess>("QstlsGuess")
      .add_property("wvg", &PyQstlsGuess::getWvg, &PyQstlsGuess::setWvg)
      .add_property("ssf", &PyQstlsGuess::getSsf, &PyQstlsGuess::setSsf)
      .add_property("adr", &PyQstlsGuess::getAdr, &PyQstlsGuess::setAdr)
      .add_property("matsubara",
                    &PyQstlsGuess::getMatsubara,
                    &PyQstlsGuess::setMatsubara);

  // Class for the free energy integrand of the VSStls scheme
  bp::class_<VSStlsInput::FreeEnergyIntegrand>("FreeEnergyIntegrand")
      .add_property("grid",
                    &PyFreeEnergyIntegrand::getGrid,
                    &PyFreeEnergyIntegrand::setGrid)
      .add_property("integrand",
                    &PyFreeEnergyIntegrand::getIntegrand,
                    &PyFreeEnergyIntegrand::setIntegrand)
      .add_property("alpha",
                    &PyFreeEnergyIntegrand::getAlpha,
                    &PyFreeEnergyIntegrand::setAlpha);

  // Class to solve the classical RPA scheme
  bp::class_<Rpa>("Rpa", bp::init<const RpaInput>())
      .def("compute", &PyRpa::compute)
      .def("rdf", &PyRpa::getRdf)
      .add_property("inputs", &PyRpa::getInput)
      .add_property("idr", &PyRpa::getIdr)
      .add_property("sdr", &PyRpa::getSdr)
      .add_property("slfc", &PyRpa::getSlfc)
      .add_property("ssf", &PyRpa::getSsf)
      .add_property("ssf_HF", &PyRpa::getSsfHF)
      .add_property("internal_energy", &PyRpa::getUInt)
      .add_property("wvg", &PyRpa::getWvg)
      .add_property("recovery", &PyRpa::getRecoveryFileName);

  // Class to solve the classical ESA scheme
  bp::class_<ESA, bp::bases<Rpa>>("ESA", bp::init<const RpaInput>())
      .def("compute", &ESA::compute);

  // Class to solve classical schemes
  bp::class_<Stls, bp::bases<Rpa>>("Stls", bp::init<const StlsInput>())
      .def("compute", &PyStls::compute)
      .add_property("inputs", &PyStls::getInput)
      .add_property("error", &PyStls::getError)
      .add_property("bf", &PyStls::getBf);

  // Class to solve the classical VS scheme
  bp::class_<VSStls, bp::bases<Rpa>>("VSStls", bp::init<const VSStlsInput>())
      .def("compute", &PyVSStls::compute)
      .add_property("inputs", &PyVSStls::getInput)
      .add_property("error", &PyVSStls::getError)
      .add_property("alpha", &PyVSStls::getAlpha)
      .add_property("free_energy_integrand", &PyVSStls::getFreeEnergyIntegrand)
      .add_property("free_energy_grid", &PyVSStls::getFreeEnergyGrid);

  // Class to solve quantum schemes
  bp::class_<Qstls, bp::bases<Stls>>("Qstls", bp::init<const QstlsInput>())
      .def("compute", &PyQstls::compute)
      .add_property("inputs", &PyQstls::getInput)
      .add_property("error", &PyQstls::getError)
      .add_property("adr", &PyQstls::getAdr);

  // Class to solve the quantum VS scheme
  bp::class_<QVSStls, bp::bases<Rpa>>("QVSStls", bp::init<const QVSStlsInput>())
      .def("compute", &PyQVSStls::compute)
      .add_property("inputs", &PyQVSStls::getInput)
      .add_property("error", &PyQVSStls::getError)
      .add_property("free_energy_integrand", &PyQVSStls::getFreeEnergyIntegrand)
      .add_property("free_energy_grid", &PyQVSStls::getFreeEnergyGrid)
      .add_property("adr", &PyQVSStls::getAdr)
      .add_property("alpha", &PyQVSStls::getAlpha);

  // MPI class
  bp::class_<PyMPI>("MPI")
      .def("rank", &PyMPI::rank)
      .def("is_root", &PyMPI::isRoot)
      .def("barrier", &PyMPI::barrier)
      .def("timer", &PyMPI::timer);

  // Post-process methods
  bp::def("compute_rdf", &PyThermo::computeRdf);
  bp::def("compute_internal_energy", &PyThermo::computeInternalEnergy);
  bp::def("compute_free_energy", &PyThermo::computeFreeEnergy);
}
