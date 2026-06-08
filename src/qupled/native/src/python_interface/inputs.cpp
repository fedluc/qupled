#include "input.hpp"
#include "python_interface/util.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace databaseUtil;
using namespace pythonUtil;

// -----------------------------------------------------------------
// Template helper functions
// -----------------------------------------------------------------

template <typename T>
py::array getAlphaGuess(T &in) {
  return toNdArray(in.getAlphaGuess());
}

template <typename T>
py::array getChemicalPotentialGuess(const T &in) {
  return toNdArray(in.getChemicalPotentialGuess());
}

template <typename T>
void setChemicalPotentialGuess(T &in, const py::list &muGuess) {
  in.setChemicalPotentialGuess(toVector(muGuess));
}

template <typename T>
void setAlphaGuess(T &in, const py::list &alphaGuess) {
  in.setAlphaGuess(toVector(alphaGuess));
}

template <typename T>
void exposeBaseInputProperties(py::class_<T> &cls) {
  cls.def_property("coupling", &T::getCoupling, &T::setCoupling)
      .def_property("degeneracy", &T::getDegeneracy, &T::setDegeneracy)
      .def_property("integral_strategy", &T::getInt2DScheme, &T::setInt2DScheme)
      .def_property("integral_error", &T::getIntError, &T::setIntError)
      .def_property("threads", &T::getNThreads, &T::setNThreads)
      .def_property("theory", &T::getTheory, &T::setTheory)
      .def_property("chemical_potential",
                    &getChemicalPotentialGuess<T>,
                    &setChemicalPotentialGuess<T>)
      .def_property("database_info", &T::getDatabaseInfo, &T::setDatabaseInfo)
      .def_property("matsubara", &T::getNMatsubara, &T::setNMatsubara)
      .def_property(
          "resolution", &T::getWaveVectorGridRes, &T::setWaveVectorGridRes)
      .def_property(
          "cutoff", &T::getWaveVectorGridCutoff, &T::setWaveVectorGridCutoff)
      .def_property(
          "frequency_cutoff", &T::getFrequencyCutoff, &T::setFrequencyCutoff)
      .def_property("dimension", &T::getDimension, &T::setDimension);
}

template <typename T>
void exposeIterativeInputProperties(py::class_<T> &cls) {
  exposeBaseInputProperties(cls);
  cls.def_property("error", &T::getErrMin, &T::setErrMin)
      .def_property("guess", &T::getGuess, &T::setGuess)
      .def_property("mixing", &T::getMixingParameter, &T::setMixingParameter)
      .def_property("iterations", &T::getNIter, &T::setNIter);
}

template <typename T>
void exposeQuantumInputProperties(py::class_<T> &cls) {
  exposeIterativeInputProperties(cls);
  cls.def_property("fixed_run_id", &T::getFixedRunId, &T::setFixedRunId);
}

template <typename T>
void exposeIetInputProperties(py::class_<T> &cls) {
  cls.def_property("mapping", &T::getMapping, &T::setMapping);
}

template <typename T>
void exposePimcInputProperties(py::class_<T> &cls) {
  cls.def_property("pimc_eta", &T::getPimcEta, &T::setPimcEta)
      .def_property("pimc_y_sec", &T::getPimcYSec, &T::setPimcYSec)
      .def_property(
          "pimc_a_cutoff", &T::getPimcACutoff, &T::setPimcACutoff);
}

template <typename T>
void exposeVSInputProperties(py::class_<T> &cls) {
  cls.def_property("error_alpha", &T::getErrMinAlpha, &T::setErrMinAlpha)
      .def_property("iterations_alpha", &T::getNIterAlpha, &T::setNIterAlpha)
      .def_property("alpha", &getAlphaGuess<T>, &setAlphaGuess<T>)
      .def_property("coupling_resolution",
                    &T::getCouplingResolution,
                    &T::setCouplingResolution)
      .def_property("degeneracy_resolution",
                    &T::getDegeneracyResolution,
                    &T::setDegeneracyResolution)
      .def_property("free_energy_integrand",
                    &T::getFreeEnergyIntegrand,
                    &T::setFreeEnergyIntegrand);
}

template <typename T>
void exposeVEInputProperties(py::class_<T> &cls) {
  cls.def_property("coupling_resolution",
                   &T::getCouplingResolution,
                   &T::setCouplingResolution)
      .def_property("free_energy_integrand",
                    &T::getFreeEnergyIntegrand,
                    &T::setFreeEnergyIntegrand);
}

void exposeInputClass(py::module_ &m) {
  auto cls = py::class_<Input>(m, "Input");
  cls.def(py::init<>());
  exposeBaseInputProperties(cls);
}

void exposeStlsInputClass(py::module_ &m) {
  auto cls = py::class_<StlsInput>(m, "StlsInput");
  cls.def(py::init<>());
  exposeIterativeInputProperties(cls);
}

namespace pythonRdfInput {
  py::array getGrid(const RecStlsInput &in) { return toNdArray(in.getRdfGrid()); }
  py::array getValues(const RecStlsInput &in) { return toNdArray(in.getRdf()); }

  void setGrid(RecStlsInput &in, const py::array_t<double> &grid) {
    in.setRdfGrid(toVector(grid));
  }

  void setValues(RecStlsInput &in, const py::array_t<double> &values) {
    in.setRdf(toVector(values));
  }
} // namespace pythonRdfInput

void exposeRecStlsInputClass(py::module_ &m) {
  auto cls = py::class_<RecStlsInput>(m, "RecStlsInput");
  cls.def(py::init<>());
  exposeBaseInputProperties(cls);
  cls.def_property("rdf_grid", pythonRdfInput::getGrid, pythonRdfInput::setGrid)
      .def_property("rdf", pythonRdfInput::getValues, pythonRdfInput::setValues);
}

void exposeStlsIetInputClass(py::module_ &m) {
  auto cls = py::class_<StlsIetInput>(m, "StlsIetInput");
  cls.def(py::init<>());
  exposeIterativeInputProperties(cls);
  exposeIetInputProperties(cls);
}

void exposeVSStlsInputClass(py::module_ &m) {
  auto cls = py::class_<VSStlsInput>(m, "VSStlsInput");
  cls.def(py::init<>());
  exposeIterativeInputProperties(cls);
  exposeVSInputProperties(cls);
}

void exposeQstlsInputClass(py::module_ &m) {
  auto cls = py::class_<QstlsInput>(m, "QstlsInput");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
}

void exposeQstlsPimcInputClass(py::module_ &m) {
  auto cls = py::class_<QstlsPimcInput>(m, "QstlsPimcInput");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
  exposePimcInputProperties(cls);
}

void exposeQstlsIetInputClass(py::module_ &m) {
  auto cls = py::class_<QstlsIetInput>(m, "QstlsIetInput");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
  exposeIetInputProperties(cls);
}

void exposeQVSStlsInputClass(py::module_ &m) {
  auto cls = py::class_<QVSStlsInput>(m, "QVSStlsInput");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
  exposeVSInputProperties(cls);
}

void exposeQVSStlsF0InputClass(py::module_ &m) {
  auto cls = py::class_<QVSStlsF0Input>(m, "QVSStlsF0Input");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
  exposePimcInputProperties(cls);
  exposeVSInputProperties(cls);
}

void exposeVEInputClass(py::module_ &m) {
  auto cls = py::class_<VEInput>(m, "VEInput");
  cls.def(py::init<>());
  exposeQuantumInputProperties(cls);
  exposeVEInputProperties(cls);
}

// -----------------------------------------------------------------
// Classes
// -----------------------------------------------------------------

namespace pythonDatabaseInfo {
  std::string getBlobStorage(const DatabaseInfo &db) { return db.blobStorage; }
  std::string getName(const DatabaseInfo &db) { return db.name; }
  std::string getRunTableName(const DatabaseInfo &db) {
    return db.runTableName;
  }
  int getRunId(const DatabaseInfo &db) { return db.runId; }

  void setBlobStorage(DatabaseInfo &db, const std::string &blobStorage) {
    db.blobStorage = blobStorage;
  }
  void setName(DatabaseInfo &db, const std::string &name) { db.name = name; }
  void setRunTableName(DatabaseInfo &db, const std::string &runTableName) {
    db.runTableName = runTableName;
  }
  void setRunId(DatabaseInfo &db, int runId) { db.runId = runId; }
} // namespace pythonDatabaseInfo

void exposeDatabaseInfoClass(py::module_ &m) {
  py::class_<DatabaseInfo>(m, "DatabaseInfo")
      .def(py::init<>())
      .def_property("blob_storage",
                    pythonDatabaseInfo::getBlobStorage,
                    pythonDatabaseInfo::setBlobStorage)
      .def_property(
          "name", pythonDatabaseInfo::getName, pythonDatabaseInfo::setName)
      .def_property(
          "run_id", pythonDatabaseInfo::getRunId, pythonDatabaseInfo::setRunId)
      .def_property("run_table_name",
                    pythonDatabaseInfo::getRunTableName,
                    pythonDatabaseInfo::setRunTableName);
}

namespace pythonGuess {
  py::array getWvg(const Guess &guess) { return toNdArray(guess.wvg); }
  py::array getSsf(const Guess &guess) { return toNdArray(guess.ssf); }
  py::array getLfc(const Guess &guess) { return toNdArray2D(guess.lfc); }

  void setWvg(Guess &guess, const py::array_t<double> &wvg) {
    guess.wvg = toVector(wvg);
  }
  void setSsf(Guess &guess, const py::array_t<double> &ssf) {
    guess.ssf = toVector(ssf);
  }
  void setLfc(Guess &guess, const py::array_t<double> &lfc) {
    if (lfc.shape(0) == 0) return;
    guess.lfc = toVector2D(lfc);
  }
} // namespace pythonGuess

void exposeGuessClass(py::module_ &m) {
  py::class_<Guess>(m, "Guess")
      .def(py::init<>())
      .def_property("wvg", pythonGuess::getWvg, pythonGuess::setWvg)
      .def_property("ssf", pythonGuess::getSsf, pythonGuess::setSsf)
      .def_property("lfc", pythonGuess::getLfc, pythonGuess::setLfc);
}

namespace pythonFreeEnergyIntegrand {
  py::array getGrid(const VSStlsInput::FreeEnergyIntegrand &fxc) {
    return toNdArray(fxc.grid);
  }
  py::array getIntegrand(const VSStlsInput::FreeEnergyIntegrand &fxc) {
    return toNdArray2D(fxc.integrand);
  }

  void setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
               const py::array_t<double> &grid) {
    fxc.grid = toVector(grid);
  }

  void setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                    const py::array &integrand) {
    fxc.integrand = toDoubleVector(integrand);
  }
} // namespace pythonFreeEnergyIntegrand

void exposeFreeEnergyIntegrandClass(py::module_ &m) {
  py::class_<VSStlsInput::FreeEnergyIntegrand>(m, "FreeEnergyIntegrand")
      .def(py::init<>())
      .def_property("grid",
                    pythonFreeEnergyIntegrand::getGrid,
                    pythonFreeEnergyIntegrand::setGrid)
      .def_property("integrand",
                    pythonFreeEnergyIntegrand::getIntegrand,
                    pythonFreeEnergyIntegrand::setIntegrand);
}

// -----------------------------------------------------------------
// Binding registration
// -----------------------------------------------------------------

namespace pythonWrappers {

  void exposeInputs(py::module_ &m) {
    exposeInputClass(m);
    exposeRecStlsInputClass(m);
    exposeStlsInputClass(m);
    exposeStlsIetInputClass(m);
    exposeVSStlsInputClass(m);
    exposeQstlsInputClass(m);
    exposeQstlsPimcInputClass(m);
    exposeQstlsIetInputClass(m);
    exposeQVSStlsInputClass(m);
    exposeQVSStlsF0InputClass(m);
    exposeVEInputClass(m);
    exposeDatabaseInfoClass(m);
    exposeGuessClass(m);
    exposeFreeEnergyIntegrandClass(m);
  }

} // namespace pythonWrappers
