#include "thermo/itcf.hpp"
#include "schemes/hf.hpp"
#include "schemes/rpa.hpp"
#include "util/numerics.hpp"
#include <cassert>
#include <cmath>

using namespace std;
using ItgParam = Integrator1D::Param;

namespace thermoUtil {

  // -----------------------------------------------------------------
  // ItcfNonInteracting class
  // -----------------------------------------------------------------

  ItcfNonInteracting::ItcfNonInteracting(const shared_ptr<const Input> in_,
                                         const double &x_,
                                         const double &mu_,
                                         const double &tau_,
                                         const double &yMin_,
                                         const double &yMax_,
                                         shared_ptr<Integrator1D> itg_,
                                         const double &idr0_)
      : in(in_),
        x(x_),
        mu(mu_),
        tau(tau_),
        yMin(yMin_),
        yMax(yMax_),
        itg(itg_),
        idr0(idr0_),
        res(numUtil::NaN) {}

  double ItcfNonInteracting::get() {
    assert(in->getDegeneracy() > 0.0);
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      shared_ptr<Integrator2D> itg2 =
          make_shared<Integrator2D>(in->getIntError());
      HFUtil::Ssf ssfCalc(
          in, x, mu, yMin, yMax, itg, vector<double>(), itg2, idr0);
      return ssfCalc.get();
    }
    compute(in->getDimension());
    return res;
  }

  void ItcfNonInteracting::compute3D() {
    auto func = [&](const double &y) -> double { return integrand(y); };
    itg->compute(func, ItgParam(yMin, yMax));
    res = itg->getSolution();
  }

  void ItcfNonInteracting::compute2D() {
    if (x > 0.0) {
      auto func = [&](const double &y) -> double { return integrand2D(y); };
      itg->compute(func, ItgParam(yMin, yMax));
      res = itg->getSolution();
    } else {
      const double Theta = in->getDegeneracy();
      res = Theta * (1.0 - exp(-1.0 / Theta));
    }
  }

  double ItcfNonInteracting::integrand(const double &y) const {
    const double Theta = in->getDegeneracy();
    const double y2 = y * y;
    if (x > 0.0) {
      const double xy = x * y;
      const double halfArg = xy / (2.0 * Theta);
      const double tauArg = xy / Theta * (tau - 0.5);
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
      const double halfArg = x * y / (2.0 * Theta);
      const double tauArg = x * y / Theta * (tau - 0.5);
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
  // ItcfHelper class
  // -----------------------------------------------------------------

  double ItcfHelper::ip() const {
    const double rs = in->getCoupling();
    if (in->getDimension() == dimensionsUtil::Dimension::D2) {
      return sqrt(2.0) * rs / x;
    } else {
      return 4.0 * numUtil::lambda * rs / (M_PI * x * x);
    }
  }

  // -----------------------------------------------------------------
  // Itcf class
  // -----------------------------------------------------------------

  Itcf::Itcf(const double &x_,
             const double &itcfHF_,
             span<const double> lfc_,
             const shared_ptr<const Input> in_,
             span<const double> idr_,
             const double &tau_)
      : x(x_),
        itcfHF(itcfHF_),
        lfc(lfc_),
        in(in_),
        idr(idr_),
        tau(tau_),
        res(numUtil::NaN) {}

  double Itcf::get() {
    assert(in->getDegeneracy() > 0.0);
    if (x == 0.0) return 0.0;
    if (in->getCoupling() == 0.0) return itcfHF;
    // For tau = 0, delegate to the Ssf calculation
    if (tau == 0.0) {
      RpaUtil::Ssf ssfCalc(x, itcfHF, lfc, in, idr);
      return ssfCalc.get();
    }
    compute(in->getDimension());
    return res;
  }

  void Itcf::compute3D() {
    const double Theta = in->getDegeneracy();
    const double suml = computeMatsubaraSummation();
    const ItcfHelper helper(in, x);
    res = itcfHF - 1.5 * helper.ip() * Theta * suml;
  }

  void Itcf::compute2D() {
    const double Theta = in->getDegeneracy();
    const double suml = computeMatsubaraSummation();
    const ItcfHelper helper(in, x);
    res = itcfHF - helper.ip() * Theta * suml;
  }

  double Itcf::computeMatsubaraSummation() const {
    const bool isStatic = lfc.size() == 1;
    const ItcfHelper helper(in, x);
    double suml = 0.0;
    for (size_t l = 0; l < idr.size(); ++l) {
      const double &idrl = idr[l];
      const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
      const double denom = 1.0 + helper.ip() * idrl * (1.0 - lfcl);
      const double f = idrl * idrl * (1.0 - lfcl) / denom;
      const double cosTerm = (l == 0) ? 1.0 : cos(2.0 * M_PI * l * tau);
      suml += (l == 0) ? f : 2.0 * f * cosTerm;
    }
    return suml;
  }

} // namespace thermoUtil
