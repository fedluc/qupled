#include "thermo_util.hpp"
#include "free_energy.hpp"
#include "internal_energy.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "rdf.hpp"
#include "dimensions_util.hpp"
#include <cassert>

using namespace std;

namespace thermoUtil {

  double computeInternalEnergy(const vector<double> &wvg,
                               const vector<double> &ssf,
                               const double &coupling,
                               const dimensionsUtil::Dimension &dim) {
    const shared_ptr<Interpolator1D> itp =
        make_shared<Interpolator1D>(wvg, ssf);
    const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>(1.0e-6);
    const InternalEnergy uInt(coupling, wvg.front(), wvg.back(), itp, itg, dim);
    return uInt.get();
  }

  double computeFreeEnergy(const vector<double> &grid,
                           const vector<double> &rsu,
                           const double &coupling) {
    return computeFreeEnergy(grid, rsu, coupling, true);
  }

  double computeFreeEnergy(const vector<double> &grid,
                           const vector<double> &rsu,
                           const double &coupling,
                           const bool normalize) {
    if (numUtil::largerThan(coupling, grid.back())) {
      MPIUtil::throwError(
          "The coupling parameter is out of range"
          " for the current grid, the free energy cannot be computed");
    }
    const shared_ptr<Interpolator1D> itp =
        make_shared<Interpolator1D>(grid, rsu);
    const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>(1.0e-6);
    const FreeEnergy freeEnergy(coupling, itp, itg, normalize);
    return freeEnergy.get();
  }


  std::vector<double> computeRdf(
      const std::vector<double> &r,
      const std::vector<double> &wvg,
      const std::vector<double> &ssf,
      const dimensionsUtil::Dimension &dim
  ) {
      assert(ssf.size() > 0 && wvg.size() > 0);
      const auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
      const int nr = r.size();
      std::vector<double> rdf(nr);
      const auto itg = std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-6);
      const auto itgf = std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-6);
      
      for (int i = 0; i < nr; ++i) {
          Rdf rdfTmp(r[i], wvg.back(), itp, itg, itgf, dim);
          rdf[i] = rdfTmp.get();
      }
      return rdf;
  }

} // namespace thermoUtil
