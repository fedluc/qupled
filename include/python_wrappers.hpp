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

// PyRpa

class PyRpa : virtual private Rpa {
public:
  PyRpa(const RpaInput& in) : Rpa(in) { ; }
  bn::ndarray getIdr4py() const;
  bn::ndarray getRdf4py(const bn::ndarray &r) const;
  bn::ndarray getSdr4py() const;
  bn::ndarray getSlfc4py() const;
  bn::ndarray getSsf4py() const;
  bn::ndarray getSsfHF4py() const;
  bn::ndarray getWvg4py() const;
  double getUInt4py() const;
  string getRecoveryFileName4py() const;
};

// PyESA

class PyESA : private ESA, public PyRpa {
public:
  PyESA(const RpaInput& in) : Rpa(in, true, false), ESA(in), PyRpa(in) { ; }
};

#endif
