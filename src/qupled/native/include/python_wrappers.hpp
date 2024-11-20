#ifndef PYTHON_WRAPPERS_HPP
#define PYTHON_WRAPPERS_HPP

#include "esa.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "qstls.hpp"
#include "qvs.hpp"
#include "rpa.hpp"
#include "stls.hpp"
#include "vsstls.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Wrapper for exposing the Input class to Python
// -----------------------------------------------------------------

class PyRpaInput {
public:

  static bn::ndarray getChemicalPotentialGuess(RpaInput &in);
  static void setChemicalPotentialGuess(RpaInput &in, const bp::list &muGuess);
};

// -----------------------------------------------------------------
// Wrapper for exposing the StlsGuess class to Python
// -----------------------------------------------------------------

class PyStlsGuess {
public:

  static bn::ndarray getWvg(const StlsInput::Guess &guess);
  static bn::ndarray getSlfc(const StlsInput::Guess &guess);
  static void setWvg(StlsInput::Guess &guess, const bn::ndarray &wvg);
  static void setSlfc(StlsInput::Guess &guess, const bn::ndarray &slfc);
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
  static bn::ndarray getAlpha(const VSStlsInput::FreeEnergyIntegrand &fxc);
  static void setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
                      const bn::ndarray &grid);
  static void setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                           const bn::ndarray &integrand);
  static void setAlpha(VSStlsInput::FreeEnergyIntegrand &fxc,
                       const bn::ndarray &alpha);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QstlsGuess class to Python
// -----------------------------------------------------------------

class PyQstlsGuess {
public:

  static bn::ndarray getWvg(const QstlsInput::Guess &guess);
  static bn::ndarray getSsf(const QstlsInput::Guess &guess);
  static bn::ndarray getAdr(const QstlsInput::Guess &guess);
  static int getMatsubara(const QstlsInput::Guess &guess);
  static void setWvg(QstlsInput::Guess &guess, const bn::ndarray &wvg);
  static void setSsf(QstlsInput::Guess &guess, const bn::ndarray &ssf);
  static void setAdr(QstlsInput::Guess &guess, const bn::ndarray &ssf);
  static void setMatsubara(QstlsInput::Guess &guess, const int matsubara);
};

// -----------------------------------------------------------------
// Wrapper for exposing the class Rpa class to Python
// -----------------------------------------------------------------

class PyRpa {
public:

  static int compute(Rpa &rpa);
  static RpaInput getInput(const Rpa &rpa);
  static bn::ndarray getIdr(const Rpa &rpa);
  static bn::ndarray getRdf(const Rpa &rpa, const bn::ndarray &r);
  static bn::ndarray getSdr(const Rpa &rpa);
  static bn::ndarray getSlfc(const Rpa &rpa);
  static bn::ndarray getSsf(const Rpa &rpa);
  static bn::ndarray getSsfHF(const Rpa &rpa);
  static bn::ndarray getWvg(const Rpa &rpa);
  static double getUInt(const Rpa &rpa);
  static std::string getRecoveryFileName(const Rpa &rpa);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Stls class to Python
// -----------------------------------------------------------------

class PyStls {
public:

  static int compute(Stls &stls);
  static StlsInput getInput(const Stls &stls);
  static double getError(const Stls &stls);
  static bn::ndarray getBf(const Stls &stls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStls class to Python
// -----------------------------------------------------------------

class PyVSStls {
public:

  static int compute(VSStls &vsstls);
  static VSStlsInput getInput(const VSStls &vsstls);
  static double getError(const VSStls &vsstls);
  static bn::ndarray getAlpha(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyIntegrand(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyGrid(const VSStls &vsstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Qstls class to Python
// -----------------------------------------------------------------

class PyQstls {
public:

  static int compute(Qstls &qstls);
  static QstlsInput getInput(const Qstls &qstls);
  static double getError(const Qstls &qstls);
  static bn::ndarray getAdr(const Qstls &qstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QVSStls class to Python
// -----------------------------------------------------------------

class PyQVSStls {
public:

  static int compute(QVSStls &qvsstls);
  static QVSStlsInput getInput(const QVSStls &qvsstls);
  static double getError(const QVSStls &qvsstls);
  static bn::ndarray getAlpha(const QVSStls &qvsstls);
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

#endif
