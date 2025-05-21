#include "python_util.hpp"
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

// -----------------------------------------------------------------
// WORK IN PROGRESS: Start
// -----------------------------------------------------------------

namespace newPythonWrappers {

  template <typename T>
  bn::ndarray getIdr(const T &scheme) {
    return pythonUtil::toNdArray2D(scheme.getIdr());
  }

  template <typename T>
  bn::ndarray getRdf(const T &scheme, const bn::ndarray &r) {
    return pythonUtil::toNdArray(scheme.getRdf(pythonUtil::toVector(r)));
  }

  template <typename T>
  bn::ndarray getSdr(const T &scheme) {
    return pythonUtil::toNdArray(scheme.getSdr());
  }

  template <typename T>
  bn::ndarray getLfc(const T &scheme) {
    return pythonUtil::toNdArray2D(scheme.getLfc());
  }

  template <typename T>
  bn::ndarray getSsf(const T &scheme) {
    return pythonUtil::toNdArray(scheme.getSsf());
  }

  template <typename T>
  bn::ndarray getWvg(const T &scheme) {
    return pythonUtil::toNdArray(scheme.getWvg());
  }

  template <typename T>
  bn::ndarray getBf(const T &scheme) {
    return pythonUtil::toNdArray(scheme.getBf());
  }

  template <typename T>
  bn::ndarray getFreeEnergyIntegrand(const T &scheme) {
    return pythonUtil::toNdArray2D(scheme.getFreeEnergyIntegrand());
  }

  template <typename T>
  bn::ndarray getFreeEnergyGrid(const T &scheme) {
    return pythonUtil::toNdArray(scheme.getFreeEnergyGrid());
  }

  template <typename Scheme>
  void exposeBaseSchemeProperties(bp::class_<Scheme> &cls) {
    cls.def("compute", &Scheme::compute)
        .def("rdf", &getRdf<Scheme>)
        .add_property("idr", &getIdr<Scheme>)
        .add_property("sdr", &getSdr<Scheme>)
        .add_property("lfc", &getLfc<Scheme>)
        .add_property("ssf", &getSsf<Scheme>)
        .add_property("uint", &Scheme::getUInt)
        .add_property("wvg", &getWvg<Scheme>);
  }

  template <typename Scheme>
  void exposeIterativeSchemeProperties(bp::class_<Scheme> &cls) {
    exposeBaseSchemeProperties(cls);
    cls.add_property("error", &Scheme::getError);
  }

  template <typename Scheme, typename Input>
  void exposeBaseSchemeClass(const std::string &className) {
    bp::class_<Scheme> cls(className.c_str(), bp::init<const Input>());
    exposeBaseSchemeProperties(cls);
  }

  template <typename Scheme, typename Input>
  void exposeIterativeSchemeClass(const std::string &className) {
    bp::class_<Scheme> cls(className.c_str(), bp::init<const Input>());
    exposeIterativeSchemeProperties(cls);
  }

  template <typename Scheme, typename Input>
  void exposeIetSchemeClass(const std::string &className) {
    bp::class_<Scheme> cls(className.c_str(), bp::init<const Input>());
    exposeIterativeSchemeProperties(cls);
    cls.add_property("bf", &getBf<Scheme>);
  }

  template <typename Scheme, typename Input>
  void exposeVSSchemeClass(const std::string &className) {
    bp::class_<Scheme> cls(className.c_str(), bp::init<const Input>());
    exposeIterativeSchemeProperties(cls);
    cls.add_property("alpha", &Scheme::getAlpha);
    cls.add_property("free_energy_integrand", &getFreeEnergyIntegrand<Scheme>);
    cls.add_property("free_energy_grid", &getFreeEnergyGrid<Scheme>);
  }

} // namespace newPythonWrappers

// -----------------------------------------------------------------
// WORK IN PROGRESS: End
// -----------------------------------------------------------------

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

  // Class for the input of the StlsIet scheme
  bp::class_<QstlsIetInput, bp::bases<IetInput, QstlsInput>>("QstlsIetInput");

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
  bp::class_<Guess>("Guess")
      .add_property("wvg", &PyGuess::getWvg, &PyGuess::setWvg)
      .add_property("ssf", &PyGuess::getSsf, &PyGuess::setSsf)
      .add_property("lfc", &PyGuess::getLfc, &PyGuess::setLfc);

  // Class for the free energy integrand of the VSStls scheme
  bp::class_<VSStlsInput::FreeEnergyIntegrand>("FreeEnergyIntegrand")
      .add_property("grid",
                    &PyFreeEnergyIntegrand::getGrid,
                    &PyFreeEnergyIntegrand::setGrid)
      .add_property("integrand",
                    &PyFreeEnergyIntegrand::getIntegrand,
                    &PyFreeEnergyIntegrand::setIntegrand);

  newPythonWrappers::exposeBaseSchemeClass<PyHF, Input>("HF");
  newPythonWrappers::exposeBaseSchemeClass<PyRpa, Input>("Rpa");
  newPythonWrappers::exposeBaseSchemeClass<PyESA, Input>("ESA");
  newPythonWrappers::exposeIterativeSchemeClass<PyStls, StlsInput>("Stls");
  newPythonWrappers::exposeIterativeSchemeClass<PyQstls, QstlsInput>("Qstls");
  newPythonWrappers::exposeIetSchemeClass<PyStlsIet, StlsIetInput>("StlsIet");
  newPythonWrappers::exposeIetSchemeClass<PyQstlsIet, QstlsIetInput>(
      "QstlsIet");
  newPythonWrappers::exposeVSSchemeClass<PyVSStls, VSStlsInput>("VSStls");
  newPythonWrappers::exposeVSSchemeClass<PyQVSStls, QVSStlsInput>("QVSStlsIet");

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
