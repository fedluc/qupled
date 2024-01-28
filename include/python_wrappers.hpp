#ifndef PYTHON_WRAPPERS_HPP
#define PYTHON_WRAPPERS_HPP

#include "rpa.hpp"
#include "esa.hpp"

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

class PyRpa  {
public:
  static bn::ndarray getIdr(const Rpa& rpa);
  static bn::ndarray getRdf(const Rpa&,
			    const bn::ndarray &r);
  static bn::ndarray getSdr(const Rpa& rpa);
  static bn::ndarray getSlfc(const Rpa& rpa);
  static bn::ndarray getSsf(const Rpa& rpa);
  static bn::ndarray getSsfHF(const Rpa& rpa);
  static bn::ndarray getWvg(const Rpa& rpa);
  static double getUInt(const Rpa& rpa);
  static string getRecoveryFileName(const Rpa& rpa);
};

// // -----------------------------------------------------------------
// // Wrapper for exposing the Stls class to Python
// // -----------------------------------------------------------------

// class PyStls : public PyRpa {
// protected:
//   Stls scheme;
// public:
//   PyStls(const StlsInput& in): PyRpa(in), scheme(in) { ; }
// }

#endif
