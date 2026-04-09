#include "thermo/thermo_util.hpp"
#include "thermo/free_energy.hpp"
#include "thermo/internal_energy.hpp"
#include "thermo/itcf.hpp"
#include "thermo/rdf.hpp"
#include "util/dimensions_util.hpp"
#include "util/mpi_util.hpp"
#include "util/numerics.hpp"
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

  std::vector<double> computeRdf(const std::vector<double> &r,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &ssf,
                                 const dimensionsUtil::Dimension &dim) {
    assert(ssf.size() > 0 && wvg.size() > 0);
    const auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
    const int nr = r.size();
    std::vector<double> rdf(nr);
    const auto itg =
        std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-6);
    const auto itgf =
        std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-6);

    for (int i = 0; i < nr; ++i) {
      Rdf rdfTmp(r[i], wvg.back(), itp, itg, itgf, dim);
      rdf[i] = rdfTmp.get();
    }
    return rdf;
  }

  Vector2D
  computeItcfNonInteractingGround(const std::shared_ptr<const Input> &in,
                                  const std::vector<double> &wvg,
                                  const std::vector<double> &tauValues) {
    auto itg2 = make_shared<Integrator2D>(in->getIntError());
    const size_t nx = wvg.size();
    const size_t ntau = tauValues.size();
    Vector2D result(nx, ntau);
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ntau; ++j) {
        ItcfNonInteractingGround itcfTmp(wvg[i], in, tauValues[j], itg2);
        result(i, j) = itcfTmp.get();
      }
    }
    return result;
  }

  Vector2D
  computeItcfNonInteractingFinite(const std::shared_ptr<const Input> &in,
                                  const std::vector<double> &wvg,
                                  const std::vector<double> &tauValues,
                                  const double mu,
                                  const Vector2D &idr) {
    using ItgType = Integrator1D::Type;
    auto itg = make_shared<Integrator1D>(ItgType::DEFAULT, in->getIntError());
    const size_t nx = wvg.size();
    const size_t ntau = tauValues.size();
    Vector2D result(nx, ntau);
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ntau; ++j) {
        ItcfNonInteracting itcfTmp(wvg[i],
                                   in,
                                   tauValues[j],
                                   mu,
                                   wvg.front(),
                                   wvg.back(),
                                   itg,
                                   idr(i, 0));
        result(i, j) = itcfTmp.get();
      }
    }
    return result;
  }

  Vector2D computeItcfNonInteracting(const std::shared_ptr<const Input> &in,
                                     const std::vector<double> &wvg,
                                     const std::vector<double> &tauValues,
                                     const double mu,
                                     const Vector2D &idr) {
    if (wvg.empty() || tauValues.empty()) { return Vector2D(); }
    if (in->getDegeneracy() == 0.0) {
      return computeItcfNonInteractingGround(in, wvg, tauValues);
    }
    return computeItcfNonInteractingFinite(in, wvg, tauValues, mu, idr);
  }

  void applyItcfGroundCorrection(const std::shared_ptr<const Input> &in,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &tauValues,
                                 const Vector2D &lfc,
                                 Vector2D &result) {
    auto itg = make_shared<Integrator1D>(in->getIntError());
    const size_t nx = wvg.size();
    const size_t ntau = tauValues.size();
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ntau; ++j) {
        ItcfGround itcfTmp(wvg[i], in, tauValues[j], result(i, j), lfc[i], itg);
        result(i, j) = itcfTmp.get();
      }
    }
  }

  void applyItcfFiniteCorrection(const std::shared_ptr<const Input> &in,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &tauValues,
                                 const Vector2D &idr,
                                 const Vector2D &lfc,
                                 Vector2D &result) {
    const size_t nx = wvg.size();
    const size_t ntau = tauValues.size();
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ntau; ++j) {
        Itcf itcfTmp(wvg[i], in, tauValues[j], result(i, j), lfc[i], idr[i]);
        result(i, j) = itcfTmp.get();
      }
    }
  }

  Vector2D computeItcf(const std::shared_ptr<const Input> &in,
                       const std::vector<double> &wvg,
                       const std::vector<double> &tauValues,
                       const double mu,
                       const Vector2D &idr,
                       const Vector2D &lfc) {
    if (wvg.empty() || tauValues.empty()) { return Vector2D(); }
    Vector2D result = computeItcfNonInteracting(in, wvg, tauValues, mu, idr);
    if (in->getTheory() == "HF") { return result; }
    if (wvg.size() != lfc.size(0)) {
      MPIUtil::throwError("Input array sizes must match");
    }
    if (in->getDegeneracy() == 0.0) {
      applyItcfGroundCorrection(in, wvg, tauValues, lfc, result);
    } else {
      applyItcfFiniteCorrection(in, wvg, tauValues, idr, lfc, result);
    }
    return result;
  }

} // namespace thermoUtil
