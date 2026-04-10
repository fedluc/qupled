#ifdef USE_MPI
#define OMPI_SKIP_MPICXX 1 // Disable OpenMPI C++ bindings
#define MPICH_SKIP_MPICXX 1 // Disable MPICH C++ bindings
#include <mpi.h>
#endif

#include "util/num_util.hpp"
#include "util/numerics.hpp"
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
