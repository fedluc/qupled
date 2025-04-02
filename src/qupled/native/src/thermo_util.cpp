#include "thermo_util.hpp"
#include "free_energy.hpp"
#include "internal_energy.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "rdf.hpp"
#include <cassert>

using namespace std;

namespace thermoUtil {

  double computeInternalEnergy(const vector<double> &wvg,
                               const vector<double> &ssf,
                               const double &coupling) {
    const shared_ptr<Interpolator1D> itp =
        make_shared<Interpolator1D>(wvg, ssf);
    const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>(1.0e-6);
    const InternalEnergy uInt(coupling, wvg.front(), wvg.back(), itp, itg);
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

  vector<double> computeRdf(const vector<double> &r,
                            const vector<double> &wvg,
                            const vector<double> &ssf) {
    assert(ssf.size() > 0 && wvg.size() > 0);
    const shared_ptr<Interpolator1D> itp =
        make_shared<Interpolator1D>(wvg, ssf);
    const int nr = r.size();
    vector<double> rdf(nr);
    const shared_ptr<Integrator1D> itg =
        make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-6);
    const shared_ptr<Integrator1D> itgf =
        make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-6);
    for (int i = 0; i < nr; ++i) {
      const Rdf rdfTmp(r[i], wvg.back(), itp, itg, itgf);
      rdf[i] = rdfTmp.get();
    }
    return rdf;
  }

} // namespace thermoUtil
