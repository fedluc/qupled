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
// Wrapper for exposing the StlsGuess class to Python
// -----------------------------------------------------------------

class PyStlsGuess {
public:

  static bn::ndarray getWvg(const StlsInput::Guess &guess);
  static bn::ndarray getSsf(const StlsInput::Guess &guess);
  static void setWvg(StlsInput::Guess &guess, const bn::ndarray &wvg);
  static void setSsf(StlsInput::Guess &guess, const bn::ndarray &lfc);
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
// Wrapper for exposing the QstlsGuess class to Python
// -----------------------------------------------------------------

class PyQstlsGuess {
public:

  static bn::ndarray getAdr(const QstlsInput::Guess &guess);
  static int getMatsubara(const QstlsInput::Guess &guess);
  static void setAdr(QstlsInput::Guess &guess, const bn::ndarray &ssf);
  static void setMatsubara(QstlsInput::Guess &guess, const int matsubara);
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
// Wrapper for exposing the Qstls class to Python
// -----------------------------------------------------------------

class PyQstls {
public:

  static bn::ndarray getAdr(const Qstls &qstls);
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
  static bn::ndarray getAdr(const QVSStls &qvsstls);
  static bn::ndarray getFreeEnergyIntegrand(const QVSStls &qvsstls);
  static bn::ndarray getFreeEnergyGrid(const QVSStls &qvsstls);
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
