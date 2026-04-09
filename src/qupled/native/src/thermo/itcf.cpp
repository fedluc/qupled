#include "thermo/itcf.hpp"
#include "schemes/hf.hpp"
#include "schemes/rpa.hpp"
#include "util/mpi_util.hpp"
#include "util/numerics.hpp"
#include <cassert>
#include <cmath>

using namespace std;
using ItgParam = Integrator1D::Param;

namespace thermoUtil {

  // -----------------------------------------------------------------
  // ItcfBase class
  // -----------------------------------------------------------------

  double ItcfBase::ip() const {
    const double rs = in->getCoupling();
    if (in->getDimension() == dimensionsUtil::Dimension::D2) {
      return sqrt(2.0) * rs / x;
    } else {
      return 4.0 * numUtil::lambda * rs / (M_PI * x * x);
    }
  }

  // -----------------------------------------------------------------
  // ItcfNonInteracting class
  // -----------------------------------------------------------------

  ItcfNonInteracting::ItcfNonInteracting(const double &x_,
                                         const shared_ptr<const Input> in_,
                                         const double &tau_,
                                         const double &mu_,
                                         const double &yMin_,
                                         const double &yMax_,
                                         shared_ptr<Integrator1D> itg_,
                                         const double &idr0_)
      : ItcfBase(x_, in_, tau_),
        mu(mu_),
        yMin(yMin_),
        yMax(yMax_),
        itg(itg_),
        idr0(idr0_),
        res(numUtil::NaN) {}

  double ItcfNonInteracting::get() {
    compute(in->getDimension());
    return res;
  }

  void ItcfNonInteracting::compute3D() {
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      shared_ptr<Integrator2D> itg2 =
          make_shared<Integrator2D>(in->getIntError());
      HFUtil::Ssf ssf(in, x, mu, yMin, yMax, itg, vector<double>(), itg2, idr0);
      res = ssf.get();
      return;
    }
    auto func = [&](const double &y) -> double { return integrand3D(y); };
    itg->compute(func, ItgParam(yMin, yMax));
    res = itg->getSolution();
  }

  void ItcfNonInteracting::compute2D() {
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      shared_ptr<Integrator2D> itg2 =
          make_shared<Integrator2D>(in->getIntError());
      HFUtil::Ssf ssf(in, x, mu, yMin, yMax, itg, vector<double>(), itg2, idr0);
      res = ssf.get();
      return;
    }
    if (x > 0.0) {
      auto func = [&](const double &y) -> double { return integrand2D(y); };
      itg->compute(func, ItgParam(yMin, yMax));
      res = itg->getSolution();
    } else {
      const double Theta = in->getDegeneracy();
      res = Theta * (1.0 - exp(-1.0 / Theta));
    }
  }

  double ItcfNonInteracting::integrand3D(const double &y) const {
    const double Theta = in->getDegeneracy();
    const double y2 = y * y;
    if (x > 0.0) {
      const double xy = x * y;
      const double halfArg = xy / (2.0 * Theta);
      const double tauArg = xy * tau - halfArg;
      const double ymx = y - x;
      const double ypx = y + x;
      const double logNum = mu - ymx * ymx / (4.0 * Theta);
      const double logDen = mu - ypx * ypx / (4.0 * Theta);
      const double logRatio = log((1.0 + exp(logNum)) / (1.0 + exp(logDen)));
      return 3.0 * Theta / 8.0 * cosh(tauArg) / sinh(halfArg) * logRatio;
    }
    return 1.5 * Theta / (1.0 + exp(y2 / Theta - mu));
  }

  double ItcfNonInteracting::integrand2D(const double &y) const {
    const double Theta = in->getDegeneracy();
    if (x > 0.0) {
      const double xy = x * y;
      const double halfArg = xy / (2.0 * Theta);
      const double tauArg = xy * tau - halfArg;
      const double ymx = y - x;
      const double ypx = y + x;
      const double eta1 = mu - ymx * ymx / (4.0 * Theta);
      const double eta2 = mu - ypx * ypx / (4.0 * Theta);
      const double fdDiff = SpecialFunctions::fermiDiracm12(eta1)
                            - SpecialFunctions::fermiDiracm12(eta2);
      return 0.5 * sqrt(Theta) * cosh(tauArg) / sinh(halfArg) * fdDiff / M_PI;
    }
    return 0.0;
  }

  // -----------------------------------------------------------------
  // ItcfNonInteractingGround class
  // -----------------------------------------------------------------

  ItcfNonInteractingGround::ItcfNonInteractingGround(
      const double &x_, const shared_ptr<const Input> in_, const double &tau_)
      : ItcfBase(x_, in_, tau_) {
    if (in->getDimension() == dimensionsUtil::Dimension::D2) {
      MPIUtil::throwError(
          "Ground state calculations of the imaginary-time "
          "correlation function in two dimensions are not implemented.");
    }
  }

  double ItcfNonInteractingGround::get() const {
    if (tau == 0.0) {
      auto itg = make_shared<Integrator1D>(Integrator1D::Type::DEFAULT,
                                           in->getIntError());
      HFUtil::SsfGround ssf(in, x, itg);
      return ssf.get();
    }
    const double xt = x * tau;
    const double xt3 = xt * xt * xt;
    const double fact = 3.0 / 16.0 / xt3;
    const double txtp1 = 2.0 * xt + 1.0;
    const double xtxm2 = xt * (x - 2.0);
    if (x == 0.0) {
      return 0.0;
    } else if (x < 2.0) {
      const double x2t2 = 2.0 * x * xt;
      return fact * (x2t2 - exp(xtxm2) * (1.0 - exp(-x2t2)) * txtp1);
    } else {
      const double em4xt = exp(-4.0 * xt);
      return fact * exp(-xtxm2) * (txtp1 * em4xt + 2.0 * xt - 1.0);
    }
  }

  // -----------------------------------------------------------------
  // Itcf class
  // -----------------------------------------------------------------

  Itcf::Itcf(const double &x_,
             const shared_ptr<const Input> in_,
             const double &tau_,
             const double &itcfHF_,
             span<const double> lfc_,
             span<const double> idr_)
      : ItcfBase(x_, in_, tau_),
        itcfHF(itcfHF_),
        lfc(lfc_),
        idr(idr_),
        res(numUtil::NaN) {}

  double Itcf::get() {
    if (x == 0.0) return 0.0;
    if (in->getCoupling() == 0.0) return itcfHF;
    compute(in->getDimension());
    return res;
  }

  void Itcf::compute3D() {
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      RpaUtil::Ssf ssf(x, itcfHF, lfc, in, idr);
      res = ssf.get();
      return;
    }
    const double Theta = in->getDegeneracy();
    const double suml = computeMatsubaraSummation();
    res = itcfHF - 1.5 * ip() * Theta * suml;
  }

  void Itcf::compute2D() {
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      RpaUtil::Ssf ssf(x, itcfHF, lfc, in, idr);
      res = ssf.get();
      return;
    }
    const double Theta = in->getDegeneracy();
    const double suml = computeMatsubaraSummation();
    res = itcfHF - ip() * Theta * suml;
  }

  double Itcf::computeMatsubaraSummation() const {
    const bool isStatic = lfc.size() == 1;
    const double &Theta = in->getDegeneracy();
    double suml = 0.0;
    for (size_t l = 0; l < idr.size(); ++l) {
      const double &idrl = idr[l];
      const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
      const double denom = 1.0 + ip() * idrl * (1.0 - lfcl);
      const double f = idrl * idrl * (1.0 - lfcl) / denom;
      const double cosTerm = (l == 0) ? 1.0 : cos(2.0 * M_PI * l * tau * Theta);
      suml += (l == 0) ? f : 2.0 * f * cosTerm;
    }
    return suml;
  }

  // -----------------------------------------------------------------
  // ItcfGround class
  // -----------------------------------------------------------------

  ItcfGround::ItcfGround(const double &x_,
                         const shared_ptr<const Input> in_,
                         const double &tau_,
                         const double &itcfHF_,
                         span<const double> lfc_,
                         shared_ptr<Integrator1D> itg_)
      : ItcfBase(x_, in_, tau_),
        itcfHF(itcfHF_),
        lfc(lfc_),
        itg(itg_) {}

  double ItcfGround::get() {
    if (x == 0.0) return 0.0;
    if (in->getCoupling() == 0.0) return itcfHF;
    if (tau == 0.0) {
      RpaUtil::SsfGround ssf(x, itcfHF, lfc, itg, in);
      return ssf.get();
    }
    const double OmegaMax = in->getFrequencyCutoff();
    auto func = [&](const double &Omega) -> double { return integrand(Omega); };
    itg->compute(func, ItgParam(0, OmegaMax));
    return 1.5 / M_PI * itg->getSolution() + itcfHF;
  }

  double ItcfGround::integrand(const double &Omega) const {
    const double idr = HFUtil::IdrGround(in, x, Omega).get();
    return (idr / (1.0 + ip() * idr * (1.0 - lfc[0])) - idr)
           * exp(-Omega * tau);
  }

} // namespace thermoUtil
