#include "python_schemes.hpp"
#include "python_util.hpp"
#include "python_wrappers.hpp"

using namespace pythonUtil;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Template functions to expose scheme properties to Python
// -----------------------------------------------------------------

template <typename T>
bn::ndarray getAlphaGuess(T &in) {
  return toNdArray(in.getAlphaGuess());
}

template <typename T>
bn::ndarray getChemicalPotentialGuess(const T &in) {
  return toNdArray(in.getChemicalPotentialGuess());
}

template <typename T>
void setChemicalPotentialGuess(Input &in, const bp::list &muGuess) {
  in.setChemicalPotentialGuess(toVector(muGuess));
}

template <typename T>
void setAlphaGuess(T &in, const bp::list &alphaGuess) {
  in.setAlphaGuess(toVector(alphaGuess));
}

template <typename Input>
void exposeBaseInputProperties(bp::class_<Input> &cls) {
  cls.add_property("coupling", &Input::getCoupling, &Input::setCoupling);
  cls.add_property("degeneracy", &Input::getDegeneracy, &Input::setDegeneracy);
  cls.add_property(
      "integral_strategy", &Input::getInt2DScheme, &Input::setInt2DScheme);
  cls.add_property("integral_error", &Input::getIntError, &Input::setIntError);
  cls.add_property("threads", &Input::getNThreads, &Input::setNThreads);
  cls.add_property("theory", &Input::getTheory, &Input::setTheory);
  cls.add_property("chemical_potential",
                   &getChemicalPotentialGuess<Input>,
                   &setChemicalPotentialGuess<Input>);
  cls.add_property(
      "database_info", &Input::getDatabaseInfo, &Input::setDatabaseInfo);
  cls.add_property("matsubara", &Input::getNMatsubara, &Input::setNMatsubara);
  cls.add_property(
      "resolution", &Input::getWaveVectorGridRes, &Input::setWaveVectorGridRes);
  cls.add_property("cutoff",
                   &Input::getWaveVectorGridCutoff,
                   &Input::setWaveVectorGridCutoff);
  cls.add_property("frequency_cutoff",
                   &Input::getFrequencyCutoff,
                   &Input::setFrequencyCutoff);
}

template <typename Input>
void exposeIterativeInputProperties(bp::class_<Input> &cls) {
  exposeBaseInputProperties(cls);
  cls.add_property("error", &Input::getErrMin, &Input::setErrMin);
  cls.add_property("guess", &Input::getGuess, &Input::setGuess);
  cls.add_property(
      "mixing", &Input::getMixingParameter, &Input::setMixingParameter);
  cls.add_property("iterations", &Input::getNIter, &Input::setNIter);
}

template <typename Input>
void exposeQuantumInputProperties(bp::class_<Input> &cls) {
  exposeIterativeInputProperties(cls);
  cls.add_property(
      "fixed_run_id", &Input::getFixedRunId, &Input::setFixedRunId);
}

template <typename Input>
void exposeIetInputProperties(bp::class_<Input> &cls) {
  cls.add_property("mapping", &Input::getMapping, &Input::setMapping);
}

template <typename Input>
void exposeVSInputProperties(bp::class_<Input> &cls) {
  cls.add_property(
      "error_alpha", &Input::getErrMinAlpha, &Input::setErrMinAlpha);
  cls.add_property(
      "iterations_alpha", &Input::getNIterAlpha, &Input::setNIterAlpha);
  cls.add_property("alpha", &getAlphaGuess<Input>, &setAlphaGuess<Input>);
  cls.add_property("coupling_resolution",
                   &Input::getCouplingResolution,
                   &Input::setCouplingResolution);
  cls.add_property("degeneracy_resolution",
                   &Input::getDegeneracyResolution,
                   &Input::setDegeneracyResolution);
  cls.add_property("free_energy_integrand",
                   &Input::getFreeEnergyIntegrand,
                   &Input::setFreeEnergyIntegrand);
}

void exposeInputClass() {
  bp::class_<Input> cls("Input");
  exposeBaseInputProperties(cls);
}

void exposeStlsInputClass() {
  bp::class_<StlsInput> cls("StlsInput");
  exposeIterativeInputProperties(cls);
}

void exposeStlsIetInputClass() {
  bp::class_<StlsIetInput> cls("StlsIetInput");
  exposeIterativeInputProperties(cls);
  exposeIetInputProperties(cls);
}

void exposeVSStlsInputClass() {
  bp::class_<VSStlsInput> cls("VSStlsInput");
  exposeIterativeInputProperties(cls);
  exposeVSInputProperties(cls);
}

void exposeQstlsInputClass() {
  bp::class_<QstlsInput> cls("QstlsInput");
  exposeQuantumInputProperties(cls);
}

void exposeQstlsIetInputClass() {
  bp::class_<QstlsIetInput> cls("QstlsIetInput");
  exposeIterativeInputProperties(cls);
  exposeIetInputProperties(cls);
}

void exposeQVSStlsInputClass() {
  bp::class_<QVSStlsInput> cls("QVSStlsInput");
  exposeIterativeInputProperties(cls);
  exposeVSInputProperties(cls);
}

// --------------------------------------------------------------------
// Namespaces used to collect functions to expose structures to python
// --------------------------------------------------------------------

namespace PyDatabaseInfo {
  std::string getName(const DatabaseInfo &db) { return db.name; }

  std::string getRunTableName(const DatabaseInfo &db) {
    return db.runTableName;
  }

  int getRunId(const DatabaseInfo &db) { return db.runId; }

  void setName(DatabaseInfo &db, const std::string &name) { db.name = name; }

  void setRunTableName(DatabaseInfo &db, const std::string &runTableName) {
    db.runTableName = runTableName;
  }

  void setRunId(DatabaseInfo &db, const int runId) { db.runId = runId; }
} // namespace PyDatabaseInfo

void exposeDatabaseInfoClass() {
  bp::class_<DatabaseInfo> cls("DatabaseInfo");
  cls.add_property("name", PyDatabaseInfo::getName, &PyDatabaseInfo::setName);
  cls.add_property(
      "run_id", PyDatabaseInfo::getRunId, &PyDatabaseInfo::setRunId);
  cls.add_property("run_table_name",
                   PyDatabaseInfo::getRunTableName,
                   &PyDatabaseInfo::setRunTableName);
}

namespace PyGuess {
  bn::ndarray getWvg(const Guess &guess) {
    return pythonUtil::toNdArray(guess.wvg);
  }

  bn::ndarray getSsf(const Guess &guess) {
    return pythonUtil::toNdArray(guess.ssf);
  }

  bn::ndarray getLfc(const Guess &guess) {
    return pythonUtil::toNdArray2D(guess.lfc);
  }

  void setWvg(Guess &guess, const bn::ndarray &wvg) {
    guess.wvg = pythonUtil::toVector(wvg);
  }

  void setSsf(Guess &guess, const bn::ndarray &ssf) {
    guess.ssf = pythonUtil::toVector(ssf);
  }

  void setLfc(Guess &guess, const bn::ndarray &lfc) {
    if (lfc.shape(0) == 0) { return; }
    guess.lfc = pythonUtil::toVector2D(lfc);
  }
} // namespace PyGuess

void exposeGuessClass() {
  bp::class_<Guess> cls("Guess");
  cls.add_property("wvg", &PyGuess::getWvg, &PyGuess::setWvg);
  cls.add_property("ssf", &PyGuess::getSsf, &PyGuess::setSsf);
  cls.add_property("lfc", &PyGuess::getLfc, &PyGuess::setLfc);
}

namespace PyFreeEnergyIntegrand {
  bn::ndarray getGrid(const VSStlsInput::FreeEnergyIntegrand &fxc) {
    return pythonUtil::toNdArray(fxc.grid);
  }

  bn::ndarray getIntegrand(const VSStlsInput::FreeEnergyIntegrand &fxc) {
    return pythonUtil::toNdArray2D(fxc.integrand);
  }

  void setGrid(VSStlsInput::FreeEnergyIntegrand &fxc, const bn::ndarray &grid) {
    fxc.grid = pythonUtil::toVector(grid);
  }

  void setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                    const bn::ndarray &integrand) {
    fxc.integrand = pythonUtil::toDoubleVector(integrand);
  }

} // namespace PyFreeEnergyIntegrand

void exposeFreeEnergyIntegrandClass() {
  bp::class_<VSStlsInput::FreeEnergyIntegrand> cls("FreeEnergyIntegrand");
  cls.add_property(
      "grid", &PyFreeEnergyIntegrand::getGrid, &PyFreeEnergyIntegrand::setGrid);
  cls.add_property("integrand",
                   &PyFreeEnergyIntegrand::getIntegrand,
                   &PyFreeEnergyIntegrand::setIntegrand);
}

// -----------------------------------------------------------------
// All inputs classes exposed to Python
// -----------------------------------------------------------------

void exposeInputs() {
  exposeInputClass();
  exposeStlsInputClass();
  exposeStlsIetInputClass();
  exposeVSStlsInputClass();
  exposeQstlsInputClass();
  exposeQstlsIetInputClass();
  exposeQVSStlsInputClass();
  exposeDatabaseInfoClass();
  exposeGuessClass();
  exposeFreeEnergyIntegrandClass();
}