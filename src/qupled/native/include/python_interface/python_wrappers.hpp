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
