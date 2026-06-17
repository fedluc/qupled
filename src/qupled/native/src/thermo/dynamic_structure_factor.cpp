#include "thermo/dynamic_structure_factor.hpp"
#include "schemes/input.hpp"
#include "thermo/dsf_integration.hpp"
#include "util/mpi_util.hpp"
#include "util/num_util.hpp"
#include <algorithm>
#include <cmath>

using namespace std;
using ItgParam = Integrator1D::Param;

namespace {

  double fermi(const double &arg) {
    if (arg >= 0.0) {
      const double e = exp(-arg);
      return e / (1.0 + e);
    }
    return 1.0 / (1.0 + exp(arg));
  }

  double softplus(const double &x) {
    if (x > 0.0) return x + log1p(exp(-x));
    return log1p(exp(x));
  }

} // namespace

namespace thermoUtil {

  IdealDynamicResponse::IdealDynamicResponse(const double &x_,
                                             const double &Omega_,
                                             const double &Theta_,
                                             const double &mu_,
                                             shared_ptr<Integrator1D> itg_)
      : x(x_),
        Omega(Omega_),
        Theta(Theta_),
        mu(mu_),
        itg(std::move(itg_)) {
    if (x <= 0.0) MPIUtil::throwError("Wave vectors must be larger than zero");
    if (Omega < 0.0) MPIUtil::throwError("Frequencies cannot be negative");
    if (Theta <= 0.0) {
      MPIUtil::throwError(
          "The ideal dynamic response requires finite temperature");
    }
  }

  double IdealDynamicResponse::real() const {
    if (Omega == 0.0) {
      auto func = [&](const double &y) { return realStaticIntegrand(y); };
      return integrateInfinite(func, {0.5 * x}) / (Theta * x);
    }

    const double a = 0.5 * x;
    const double b = 0.5 * Omega / x;
    vector<double> singularities = {abs(a - b), a + b};
    if (b > a) singularities.push_back(b - a);
    auto func = [&](const double &y) { return realDynamicIntegrand(y); };
    return -0.5 * integrateInfinite(func, singularities) / x;
  }

  double IdealDynamicResponse::imaginary() const {
    if (Omega == 0.0) return 0.0;
    const double a = 0.5 * x;
    const double b = 0.5 * Omega / x;
    const double lower = (a - b) * (a - b);
    const double upper = (a + b) * (a + b);
    const double integral =
        Theta * (softplus(mu - lower / Theta) - softplus(mu - upper / Theta));
    return M_PI * integral / (4.0 * x);
  }

  double IdealDynamicResponse::realDynamicIntegrand(const double &y) const {
    if (y == 0.0) return 0.0;
    const double a = 0.5 * x;
    const double b = 0.5 * Omega / x;
    const double numerator = (y - a) * (y - a) - b * b;
    const double denominator = (y + a) * (y + a) - b * b;
    const double logRatio = log(abs(numerator)) - log(abs(denominator));
    return y * fermi(y * y / Theta - mu) * logRatio;
  }

  double IdealDynamicResponse::realStaticIntegrand(const double &y) const {
    const double y2 = y * y;
    const double xHalf = 0.5 * x;
    const double coeff = y2 - xHalf * xHalf;
    double logarithmicTerm = 0.0;
    if (coeff != 0.0) {
      logarithmicTerm = coeff * (log(2.0 * y + x) - log(abs(2.0 * y - x)));
    }
    const double occupation = fermi(y2 / Theta - mu);
    return (logarithmicTerm + x * y) * y * occupation * (1.0 - occupation);
  }

  double IdealDynamicResponse::integrateInfinite(
      const function<double(double)> &func,
      const vector<double> &singularities) const {
    vector<double> splitPoints = {0.0, 1.0};
    for (const double y : singularities) {
      if (y > 0.0 && isfinite(y)) splitPoints.push_back(y / (1.0 + y));
    }
    sort(splitPoints.begin(), splitPoints.end());
    splitPoints.erase(unique(splitPoints.begin(), splitPoints.end()),
                      splitPoints.end());

    auto transformed = [&](const double &t) {
      if (t <= 0.0) return func(0.0);
      if (t >= 1.0) return 0.0;
      const double inv = 1.0 / (1.0 - t);
      return func(t * inv) * inv * inv;
    };

    double result = 0.0;
    for (size_t i = 1; i < splitPoints.size(); ++i) {
      itg->compute(transformed, ItgParam(splitPoints[i - 1], splitPoints[i]));
      result += itg->getSolution();
    }
    return result;
  }

  pair<Vector2D, Vector2D>
  computeIdealDynamicResponse(const vector<double> &wvg,
                              const vector<double> &frequency,
                              const double &Theta,
                              const double &mu,
                              const double &intError) {
    if (intError <= 0.0) {
      MPIUtil::throwError("The integration accuracy must be larger than zero");
    }
    Vector2D realPart(wvg.size(), frequency.size());
    Vector2D imaginaryPart(wvg.size(), frequency.size());
    auto itg = make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, intError);
    for (size_t i = 0; i < wvg.size(); ++i) {
      for (size_t j = 0; j < frequency.size(); ++j) {
        IdealDynamicResponse response(wvg[i], frequency[j], Theta, mu, itg);
        realPart(i, j) = response.real();
        imaginaryPart(i, j) = response.imaginary();
      }
    }
    return {realPart, imaginaryPart};
  }

  double computeDSFValue(const double &x,
                         const double &Omega,
                         const double &Theta,
                         const double &mu,
                         const double &interaction,
                         const shared_ptr<Integrator1D> &itg) {
    IdealDynamicResponse response(x, Omega, Theta, mu, itg);
    const double realPart = response.real();
    const double reDenom = 1.0 + interaction * realPart;
    if (Omega == 0.0) {
      const double occupation = fermi(x * x / (4.0 * Theta) - mu);
      return Theta * occupation / (4.0 * x * reDenom * reDenom);
    }
    const double imaginaryPart = response.imaginary();
    const double imDenom = interaction * imaginaryPart;
    const double dielectric = reDenom * reDenom + imDenom * imDenom;
    return imaginaryPart / (M_PI * -expm1(-Omega / Theta) * dielectric);
  }

  Vector2D computeDSF(const shared_ptr<const Input> &in,
                      const vector<double> &wvg,
                      const vector<double> &frequency,
                      const double &mu,
                      const Vector2D &lfc) {
    if (in->getDimension() != dimensionsUtil::Dimension::D3) {
      MPIUtil::throwError(
          "The static dynamic structure factor is only implemented in 3D");
    }
    if (in->getDegeneracy() <= 0.0) {
      MPIUtil::throwError(
          "The static dynamic structure factor requires finite temperature");
    }
    if (lfc.size(0) != wvg.size() || lfc.size(1) != 1) {
      MPIUtil::throwError(
          "The dynamic structure factor requires a static local-field "
          "correction with shape (nx, 1)");
    }
    for (const double Omega : frequency) {
      if (Omega < 0.0) MPIUtil::throwError("Frequencies cannot be negative");
    }

    const double Theta = in->getDegeneracy();
    const double rs = in->getCoupling();
    Vector2D result(wvg.size(), frequency.size());
    auto itg = make_shared<Integrator1D>(Integrator1D::Type::DEFAULT,
                                         in->getIntError());
    for (size_t i = 0; i < wvg.size(); ++i) {
      const double x = wvg[i];
      if (x == 0.0) {
        result.fill(i, 0.0);
        continue;
      }
      const double interaction =
          4.0 * numUtil::lambda * rs * (1.0 - lfc(i, 0)) / (M_PI * x * x);
      for (size_t j = 0; j < frequency.size(); ++j) {
        result(i, j) = computeDSFValue(
            x, frequency[j], Theta, mu, interaction, itg);
      }
    }
    return result;
  }

  Vector2D computeAdaptiveITCFfromDSF(
      const shared_ptr<const Input> &in,
      const vector<double> &wvg,
      const vector<double> &frequency,
      const vector<double> &tau,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf) {
    const double Theta = in->getDegeneracy();
    if (Theta <= 0.0) {
      MPIUtil::throwError(
          "Adaptive DSF integration requires finite temperature");
    }
    const double beta = 1.0 / Theta;
    vector<function<double(double)>> kernels;
    for (const double t : tau) {
      if (t < 0.0 || t > beta) {
        MPIUtil::throwError("Imaginary time must lie in [0, 1/Theta]");
      }
      kernels.emplace_back([=](const double &Omega) {
        return exp(-Omega * t) + exp(-Omega * (beta - t));
      });
    }
    return computeAdaptiveWeightedDSFIntegral(
        in, wvg, frequency, mu, lfc, dsf, kernels, 4.0);
  }

} // namespace thermoUtil
