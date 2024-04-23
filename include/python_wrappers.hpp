#ifndef PYTHON_WRAPPERS_HPP
#define PYTHON_WRAPPERS_HPP

#include "rpa.hpp"
#include "esa.hpp"
#include "stls.hpp"
#include "vsstls.hpp"
#include "qstls.hpp"
#include "qvs.hpp"

// Forward declarations
namespace boost {
  namespace python {
    namespace numpy {
      class ndarray;
    }
  }
}
namespace bp = boost::python;
namespace bn = boost::python::numpy;

// -----------------------------------------------------------------
// Wrapper for exposing the Input class to Python
// -----------------------------------------------------------------

class PyRpaInput {
public:
  static bn::ndarray getChemicalPotentialGuess(RpaInput &in);
  static void setChemicalPotentialGuess(RpaInput &in,
					const bp::list &muGuess);
};

// -----------------------------------------------------------------
// Wrapper for exposing the SlfcGuess class to Python
// -----------------------------------------------------------------

class PySlfcGuess {
public:
  static bn::ndarray getWvg(const StlsInput::SlfcGuess &guess);
  static bn::ndarray getSlfc(const StlsInput::SlfcGuess &guess);
  static void setWvg(StlsInput::SlfcGuess &guess,
		     const bn::ndarray& wvg);
  static void setSlfc(StlsInput::SlfcGuess &guess,
		      const bn::ndarray& slfc);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStlsInput class to Python
// -----------------------------------------------------------------

class PyVSStlsInput {
public:
  static bn::ndarray getAlphaGuess(VSStlsInput &in);
  static void setAlphaGuess(VSStlsInput &in,
			    const bp::list &alphaGuess);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QVSStlsInput class to Python
// -----------------------------------------------------------------

class PyQVSStlsInput {
public:
  static bn::ndarray getAlphaGuess(QVSStlsInput &in);
  static void setAlphaGuess(QVSStlsInput &in,
			    const bp::list &alphaGuess);
};


// -----------------------------------------------------------------
// Wrapper for exposing the FreeEnergyIntegrand class to Python
// -----------------------------------------------------------------

class PyFreeEnergyIntegrand {
public:
  static bn::ndarray getGrid(const VSStlsInput::FreeEnergyIntegrand& fxc);
  static bn::ndarray getIntegrand(const VSStlsInput::FreeEnergyIntegrand& fxc);
  static void setGrid(VSStlsInput::FreeEnergyIntegrand& fxc,
		      const bn::ndarray& grid);
  static void setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
			   const bn::ndarray& integrand);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QstlsGuess class to Python
// -----------------------------------------------------------------

class PyQstlsGuess {
public:
  static bn::ndarray getWvg(const QstlsInput::QstlsGuess &guess);
  static bn::ndarray getSsf(const QstlsInput::QstlsGuess &guess);
  static bn::ndarray getAdr(const QstlsInput::QstlsGuess &guess);
  static int getMatsubara(const QstlsInput::QstlsGuess &guess);
  static void setWvg(QstlsInput::QstlsGuess &guess,
		     const bn::ndarray& wvg);
  static void setSsf(QstlsInput::QstlsGuess &guess,
		      const bn::ndarray& ssf);
  static void setAdr(QstlsInput::QstlsGuess &guess,
		     const bn::ndarray& ssf);
  static void setMatsubara(QstlsInput::QstlsGuess &guess,
			   const int matsubara);
};

// -----------------------------------------------------------------
// Wrapper for exposing the class Rpa class to Python
// -----------------------------------------------------------------

class PyRpa {
public:
  static int compute(Rpa& rpa);
  static bn::ndarray getIdr(const Rpa& rpa);
  static bn::ndarray getRdf(const Rpa& rpa,
			    const bn::ndarray &r);
  static bn::ndarray getSdr(const Rpa& rpa);
  static bn::ndarray getSlfc(const Rpa& rpa);
  static bn::ndarray getSsf(const Rpa& rpa);
  static bn::ndarray getSsfHF(const Rpa& rpa);
  static bn::ndarray getWvg(const Rpa& rpa);
  static double getUInt(const Rpa& rpa);
  static std::string getRecoveryFileName(const Rpa& rpa);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Stls class to Python
// -----------------------------------------------------------------

class PyStls {
public:
  static int compute(Stls& stls);
  static double getError(const Stls& stls);
  static bn::ndarray getBf(const Stls& stls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStls class to Python
// -----------------------------------------------------------------

class PyVSStls {
public:
  static int compute(VSStls& vsstls);
  static double getError(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyIntegrand(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyGrid(const VSStls &vsstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Qstls class to Python
// -----------------------------------------------------------------

class PyQstls {
public:
  static int compute(Qstls& qstls);
  static double getError(const Qstls& qstls);
  static bn::ndarray getAdr(const Qstls& qstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the QVSStls class to Python
// -----------------------------------------------------------------

class PyQVSStls {
public:
  static int compute(QVSStls& qvsstls);
  static double getError(const QVSStls &qvsstls);
  static bn::ndarray getAdr(const QVSStls& qvsstls);
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
  static int rank() { return parallelUtil::MPI::rank(); }
  static bool isRoot() { return parallelUtil::MPI::isRoot(); }
  static void barrier() { return parallelUtil::MPI::barrier(); }
  static double timer() { return parallelUtil::MPI::timer(); }
};

#endif
