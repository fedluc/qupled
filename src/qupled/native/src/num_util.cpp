#ifdef USE_MPI
#include <mpi.h>
#define OMPI_SKIP_MPICXX 1 // Disable MPI-C++ bindings
#endif

#include "num_util.hpp"
#include "numerics.hpp"
#include <cassert>
#include <numeric>
#include <omp.h>

namespace numUtil {

  bool isZero(const double &x) { return abs(x) < dtol; }

  bool equalTol(const double &x, const double &y) {
    return abs(x - y) < x * dtol;
  }

  bool largerThan(const double &x, const double &y) { return x - y > x * dtol; }

} // namespace numUtil
