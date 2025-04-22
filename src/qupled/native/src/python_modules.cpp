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
  bp::class_<Input>("Input")
      .add_property("coupling", &Input::getCoupling, &Input::setCoupling)
      .add_property("degeneracy", &Input::getDegeneracy, &Input::setDegeneracy)
      .add_property(
          "integral_strategy", &Input::getInt2DScheme, &Input::setInt2DScheme)
      .add_property("integral_error", &Input::getIntError, &Input::setIntError)
      .add_property("threads", &Input::getNThreads, &Input::setNThreads)
      .add_property("theory", &Input::getTheory, &Input::setTheory)
      .add_property("chemical_potential",
                    &PyInput::getChemicalPotentialGuess,
                    &PyInput::setChemicalPotentialGuess)
      .add_property(
          "database_info", &Input::getDatabaseInfo, &Input::setDatabaseInfo)
      .add_property("matsubara", &Input::getNMatsubara, &Input::setNMatsubara)
      .add_property("resolution",
                    &Input::getWaveVectorGridRes,
                    &Input::setWaveVectorGridRes)
      .add_property("cutoff",
                    &Input::getWaveVectorGridCutoff,
                    &Input::setWaveVectorGridCutoff)
      .add_property("frequency_cutoff",
                    &Input::getFrequencyCutoff,
                    &Input::setFrequencyCutoff);

  // Class for the input of the Stls scheme
  bp::class_<StlsInput, bp::bases<Input>>("StlsInput")
      .add_property("error", &StlsInput::getErrMin, &StlsInput::setErrMin)
      .add_property("guess", &StlsInput::getGuess, &StlsInput::setGuess)
      .add_property("mixing",
                    &StlsInput::getMixingParameter,
                    &StlsInput::setMixingParameter)
      .add_property("iterations", &StlsInput::getNIter, &StlsInput::setNIter);

  // Class for the input of the IET schemes
  bp::class_<IetInput>("IetInput")
      .add_property("mapping", &IetInput::getMapping, &IetInput::setMapping);

  // Class for the input of the StlsIet scheme
  bp::class_<StlsIetInput, bp::bases<IetInput, StlsInput>>("StlsIetInput");

  // Class for the input of the VS scheme
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
                    &VSInput::setFreeEnergyIntegrand);

  // Class for the input of the VSStls scheme
  bp::class_<VSStlsInput, bp::bases<VSInput, StlsInput>>("VSStlsInput");

  // Class for the input of the Qstls scheme
  bp::class_<QstlsInput, bp::bases<StlsInput>>("QstlsInput")
      .add_property("guess", &QstlsInput::getGuess, &QstlsInput::setGuess)
      .add_property("fixed_run_id",
                    &QstlsInput::getFixedRunId,
                    &QstlsInput::setFixedRunId)
      .add_property(
          "fixed_iet", &QstlsInput::getFixedIet, &QstlsInput::setFixedIet);

  // Class for the input of the QVSStls scheme
  bp::class_<QVSStlsInput, bp::bases<VSInput, QstlsInput>>("QVSStlsInput");

  // Class for the database information
  bp::class_<DatabaseInfo>("DatabaseInfo")
      .add_property("name", PyDatabaseInfo::getName, &PyDatabaseInfo::setName)
      .add_property(
          "run_id", PyDatabaseInfo::getRunId, &PyDatabaseInfo::setRunId)
      .add_property("run_table_name",
                    PyDatabaseInfo::getRunTableName,
                    &PyDatabaseInfo::setRunTableName);

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
                    &PyFreeEnergyIntegrand::setIntegrand);

  // Class to solve the classical RPA scheme
  bp::class_<HF>("HF", bp::init<const Input>())
      .def("compute", &HF::compute)
      .def("rdf", &pyHF::getRdf)
      .add_property("inputs", &pyHF::getInput)
      .add_property("idr", &pyHF::getIdr)
      .add_property("sdr", &pyHF::getSdr)
      .add_property("slfc", &pyHF::getSlfc)
      .add_property("ssf", &pyHF::getSsf)
      .add_property("uint", &pyHF::getUInt)
      .add_property("wvg", &pyHF::getWvg);

  // Class to solve the classical RPA scheme
  bp::class_<Rpa, bp::bases<HF>>("Rpa", bp::init<const Input>())
      .def("compute", &Rpa::compute);

  // Class to solve the classical ESA scheme
  bp::class_<ESA, bp::bases<Rpa>>("ESA", bp::init<const Input>())
      .def("compute", &ESA::compute);

  // Class to solve classical schemes
  bp::class_<Stls, bp::bases<Rpa>>("Stls", bp::init<const StlsInput>())
      .def("compute", &PyStls::compute)
      .add_property("inputs", &PyStls::getInput)
      .add_property("error", &PyStls::getError);

  // Class to solve classical schemes
  bp::class_<StlsIet, bp::bases<Stls>>("StlsIet",
                                       bp::init<const StlsIetInput>())
      .def("compute", &PyStlsIet::compute)
      .add_property("bf", &PyStlsIet::getBf);

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
      .staticmethod("rank")
      .def("is_root", &PyMPI::isRoot)
      .staticmethod("is_root")
      .def("barrier", &PyMPI::barrier)
      .staticmethod("barrier")
      .def("timer", &PyMPI::timer)
      .staticmethod("timer");

  // Post-process methods
  bp::def("compute_rdf", &PyThermo::computeRdf);
  bp::def("compute_internal_energy", &PyThermo::computeInternalEnergy);
  bp::def("compute_free_energy", &PyThermo::computeFreeEnergy);
}
