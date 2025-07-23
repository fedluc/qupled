#include "dimensions_util.hpp"
#include "mpi_util.hpp"

namespace dimensionsUtil {

  void DimensionsHandler::compute(const Dimension &dim) {
    switch (dim) {
    case Dimension::D2: compute2D(); break;
    case Dimension::D3: compute3D(); break;
    default: MPIUtil::throwError("Unsupported dimension");
    }
  }

} // namespace dimensionsUtil