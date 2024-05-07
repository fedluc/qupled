#ifdef USE_MPI
#include <mpi.h>
#define OMPI_SKIP_MPICXX 1 // Disable MPI-C++ bindings
#endif

#include <cassert>
#include <numeric>
#include <omp.h>
#include "numerics.hpp"
#include "num_util.hpp"

namespace numUtil {

  bool equalTol(const double &x, const double &y) {
    return abs(x - y) < x * dtol;
  }

  bool largerThan(const double &x, const double &y) { return x - y > x * dtol; }

} // namespace numUtil
