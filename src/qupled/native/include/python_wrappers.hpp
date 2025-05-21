#ifndef PYTHON_WRAPPERS_HPP
#define PYTHON_WRAPPERS_HPP

#include "database.hpp"
#include "esa.hpp"
#include "hf.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "qstls.hpp"
#include "qstlsiet.hpp"
#include "qvsstls.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "stlsiet.hpp"
#include "vsstls.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Wrapper for exposing the Input class to Python
// -----------------------------------------------------------------

class PyInput {
public:

  static bn::ndarray getChemicalPotentialGuess(Input &in);
  static void setChemicalPotentialGuess(Input &in, const bp::list &muGuess);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Guess class to Python
// -----------------------------------------------------------------

class PyGuess {
public:

  static bn::ndarray getWvg(const Guess &guess);
  static bn::ndarray getSsf(const Guess &guess);
  static bn::ndarray getLfc(const Guess &guess);
  static void setWvg(Guess &guess, const bn::ndarray &wvg);
  static void setSsf(Guess &guess, const bn::ndarray &lfc);
  static void setLfc(Guess &guess, const bn::ndarray &lfc);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStlsInput class to Python
// -----------------------------------------------------------------

class PyVSInput {
public:

  static bn::ndarray getAlphaGuess(VSInput &in);
  static void setAlphaGuess(VSInput &in, const bp::list &alphaGuess);
};

// -----------------------------------------------------------------
// Wrapper for exposing the FreeEnergyIntegrand class to Python
// -----------------------------------------------------------------

class PyFreeEnergyIntegrand {
public:

  static bn::ndarray getGrid(const VSStlsInput::FreeEnergyIntegrand &fxc);
  static bn::ndarray getIntegrand(const VSStlsInput::FreeEnergyIntegrand &fxc);
  static void setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
                      const bn::ndarray &grid);
  static void setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                           const bn::ndarray &integrand);
};

// -----------------------------------------------------------------
// Wrapper for exposing the class Rpa class to Python
// -----------------------------------------------------------------

class pyHF {
public:

  static bn::ndarray getIdr(const HF &hf);
  static bn::ndarray getRdf(const HF &hf, const bn::ndarray &r);
  static bn::ndarray getSdr(const HF &hf);
  static bn::ndarray getLfc(const HF &hf);
  static bn::ndarray getSsf(const HF &hf);
  static bn::ndarray getWvg(const HF &hf);
  static double getUInt(const HF &hf);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Stls class to Python
// -----------------------------------------------------------------

class PyStls {
public:

  static double getError(const Stls &stls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the StlsIet class to Python
// -----------------------------------------------------------------

class PyStlsIet {
public:

  static bn::ndarray getBf(const StlsIet &stlsiet);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStls class to Python
// -----------------------------------------------------------------

class PyVSStls {
public:

  static int compute(VSStls &vsstls);
  static double getError(const VSStls &vsstls);
  static double getAlpha(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyIntegrand(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyGrid(const VSStls &vsstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QstlsIet class to Python
// -----------------------------------------------------------------

class PyQstlsIet {
public:

  static bn::ndarray getBf(const QstlsIet &qstlsiet);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QVSStls class to Python
// -----------------------------------------------------------------

class PyQVSStls {
public:

  static int compute(QVSStls &qvsstls);
  static double getError(const QVSStls &qvsstls);
  static double getAlpha(const QVSStls &qvsstls);
  static bn::ndarray getFreeEnergyIntegrand(const QVSStls &qvsstls);
  static bn::ndarray getFreeEnergyGrid(const QVSStls &qvsstls);
};

// -----------------------------------------------------------------
// New wrappers
// -----------------------------------------------------------------

class pyHFNew : public HF {
public:

  explicit pyHFNew(const Input &in)
      : HF(std::make_shared<Input>(in)) {}
};
class pyRpaNew : public Rpa {
public:

  explicit pyRpaNew(const Input &in)
      : Rpa(std::make_shared<Input>(in)) {}
};

class pyESANew : public ESA {
public:

  explicit pyESANew(const Input &in)
      : ESA(std::make_shared<Input>(in)) {}
};

class pyStlsNew : public Stls {
public:

  explicit pyStlsNew(const StlsInput &in)
      : Stls(std::make_shared<StlsInput>(in)) {}
};

class pyStlsIetNew : public StlsIet {
public:

  explicit pyStlsIetNew(const StlsIetInput &in)
      : StlsIet(std::make_shared<StlsIetInput>(in)) {}
};

class pyVSStlsNew : public VSStls {
public:

  explicit pyVSStlsNew(const VSStlsInput &in)
      : VSStls(std::make_shared<VSStlsInput>(in)) {}
};
class pyQstlsNew : public Qstls {
public:

  explicit pyQstlsNew(const QstlsInput &in)
      : Qstls(std::make_shared<QstlsInput>(in)) {}
};

class pyQstlsIetNew : public QstlsIet {
public:

  explicit pyQstlsIetNew(const QstlsIetInput &in)
      : QstlsIet(std::make_shared<QstlsIetInput>(in)) {}
};

class pyQVSStlsNew : public QVSStls {
public:

  explicit pyQVSStlsNew(const QVSStlsInput &in)
      : QVSStls(std::make_shared<QVSStlsInput>(in)) {}
};

// -----------------------------------------------------------------
// Wrapper for exposing methods in thermoUtil to Python
// -----------------------------------------------------------------

class PyThermo {
public:

  static bn::ndarray computeRdf(const bn::ndarray &rIn,
                                const bn::ndarray &wvgIn,
                                const bn::ndarray &ssfIn);
  static double computeInternalEnergy(const bn::ndarray &wvgIn,
                                      const bn::ndarray &ssfIn,
                                      const double &coupling);
  static double computeFreeEnergy(const bn::ndarray &gridIn,
                                  const bn::ndarray &rsuIn,
                                  const double &coupling);
};
// -----------------------------------------------------------------
// Wrapper for exposing MPI methods to Python
// -----------------------------------------------------------------

class PyMPI {
public:

  static int rank() { return MPIUtil::rank(); }
  static bool isRoot() { return MPIUtil::isRoot(); }
  static void barrier() { return MPIUtil::barrier(); }
  static double timer() { return MPIUtil::timer(); }
};

// -----------------------------------------------------------------
// Wrapper for exposing MPI methods to Python
// -----------------------------------------------------------------

class PyDatabaseInfo {
public:

  static std::string getName(const DatabaseInfo &dbInfo);
  static int getRunId(const DatabaseInfo &dbInfo);
  static std::string getRunTableName(const DatabaseInfo &dbInfo);
  static void setName(DatabaseInfo &dbInfo, const std::string &name);
  static void setRunId(DatabaseInfo &dbInfo, const int runId);
  static void setRunTableName(DatabaseInfo &dbInfo,
                              const std::string &runTableName);
};

#endif
