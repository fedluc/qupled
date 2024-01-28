#ifndef PYTHON_WRAPPERS_HPP
#define PYTHON_WRAPPERS_HPP

#include "rpa.hpp"
#include "esa.hpp"
#include "stls.hpp"
#include "vsstls.hpp"
#include "qstls.hpp"

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
// Wrapper for exposing the class Rpa class to Python
// -----------------------------------------------------------------

class PyRpa {
public:
  static bn::ndarray getIdr(const Rpa& rpa);
  static bn::ndarray getRdf(const Rpa& rpa,
			    const bn::ndarray &r);
  static bn::ndarray getSdr(const Rpa& rpa);
  static bn::ndarray getSlfc(const Rpa& rpa);
  static bn::ndarray getSsf(const Rpa& rpa);
  static bn::ndarray getSsfHF(const Rpa& rpa);
  static bn::ndarray getWvg(const Rpa& rpa);
  static double getUInt(const Rpa& rpa);
  static string getRecoveryFileName(const Rpa& rpa);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Stls class to Python
// -----------------------------------------------------------------

class PyStls {
public:
  static int compute(Stls& stls);
  static bn::ndarray getBf(const Stls& stls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the VSStls class to Python
// -----------------------------------------------------------------

class PyVSStls {
public:
  static int compute(VSStls& vsstls);
  static bn::ndarray getFreeEnergyIntegrand(const VSStls &vsstls);
  static bn::ndarray getFreeEnergyGrid(const VSStls &vsstls);
};

// -----------------------------------------------------------------
// Wrapper for exposing the Qstls class to Python
// -----------------------------------------------------------------

class PyQstls {
public:
  static int compute(Qstls& qstls);
  static bn::ndarray getAdr(const Qstls& qstls);
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

#endif
