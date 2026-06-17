#include "thermo/dsf_integration.hpp"
#include "schemes/input.hpp"
#include "thermo/dynamic_structure_factor.hpp"
#include "util/mpi_util.hpp"
#include "util/num_util.hpp"
#include "util/numerics.hpp"
#include <algorithm>
#include <cmath>

using namespace std;
using ItgParam = Integrator1D::Param;

namespace thermoUtil {

  Vector2D computeAdaptiveWeightedDSFIntegral(
      const shared_ptr<const Input> &in,
      const vector<double> &wvg,
      const vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf,
      const vector<function<double(double)>> &kernels,
      const double &maxWidthFactor,
      const DSFIntegrationOptions &options) {
    const double Theta = in->getDegeneracy();
    if (Theta <= 0.0) {
      MPIUtil::throwError(
          "Adaptive DSF integration requires finite temperature");
    }
    if (frequency.size() < 3 || frequency.front() != 0.0) {
      MPIUtil::throwError("The frequency grid must start at zero and contain "
                          "at least three points");
    }
    if (dsf.size(0) != wvg.size() || lfc.size(0) != wvg.size()
        || dsf.size(1) != frequency.size()) {
      MPIUtil::throwError("Invalid grid, DSF, or local-field-correction shape");
    }
    for (size_t j = 1; j < frequency.size(); ++j) {
      if (frequency[j] <= frequency[j - 1]) {
        MPIUtil::throwError("The frequency grid must be strictly increasing");
      }
    }

    const double rs = in->getCoupling();
    const double maxFrequency = frequency.back();
    Integrator1D baselineItg(Integrator1D::Type::DEFAULT, in->getIntError());
    auto responseItg = make_shared<Integrator1D>(Integrator1D::Type::DEFAULT,
                                                 in->getIntError());
    Integrator1D peakItg(Integrator1D::Type::DEFAULT, in->getIntError());
    Vector2D result(wvg.size(), kernels.size());

    for (size_t i = 0; i < wvg.size(); ++i) {
      const span<const double> spectrum = dsf[i];
      const Interpolator1D spectrumItp(
          frequency[0], spectrum[0], frequency.size());
      auto integrateBaseline = [&](const function<double(double)> &kernel,
                                   const double lower,
                                   const double upper) {
        auto integrand = [&](const double &Omega) {
          return 1.5 * spectrumItp.eval(Omega) * kernel(Omega);
        };
        baselineItg.compute(integrand, ItgParam(lower, upper));
        return baselineItg.getSolution();
      };
      for (size_t k = 0; k < kernels.size(); ++k) {
        result(i, k) =
            integrateBaseline(kernels[k], frequency.front(), maxFrequency);
      }

      const double x = wvg[i];
      if (x == 0.0 || x > options.peakSearchXMax) continue;
      const double interaction =
          4.0 * numUtil::lambda * rs * (1.0 - lfc(i, 0)) / (M_PI * x * x);
      auto denominator = [&](const double &Omega) {
        IdealDynamicResponse response(x, Omega, Theta, mu, responseItg);
        return 1.0 + interaction * response.real();
      };
      const double scanLimit =
          min(maxFrequency, x * x + options.peakSearchOmegaBuffer);
      const double scanStep =
          min(options.peakSearchStepMax,
              scanLimit / options.peakSearchMinIntervals);
      double left = 0.0;
      double leftValue = denominator(left);
      for (double right = scanStep; right <= scanLimit + 0.5 * scanStep;
           right += scanStep) {
        right = min(right, scanLimit);
        const double rightValue = denominator(right);
        if (leftValue * rightValue < 0.0) {
          BrentRootSolver rootSolver;
          rootSolver.solve(denominator, {left, right});
          const double root = rootSolver.getSolution();
          const auto upper =
              upper_bound(frequency.begin(), frequency.end(), root);
          if (upper == frequency.begin() || upper == frequency.end()) {
            left = right;
            leftValue = rightValue;
            continue;
          }
          const size_t j = distance(frequency.begin(), upper);
          const double spacing = frequency[j] - frequency[j - 1];
          const double h =
              max(options.slopeStepAbs, options.slopeStepRel * root);
          const double slope =
              (denominator(root + h) - denominator(max(0.0, root - h)))
              / (root + h - max(0.0, root - h));
          IdealDynamicResponse response(x, root, Theta, mu, responseItg);
          const double width = abs(interaction * response.imaginary() / slope);
          if (width < maxWidthFactor * spacing) {
            const double detailedBalance = -expm1(-root / Theta);
            const double peakArea =
                1.0 / (detailedBalance * interaction * abs(slope));
            const double halfWindow =
                max(options.peakWindowSpacingFactor * spacing,
                    options.peakWindowWidthFactor * width);
            const auto windowLower = lower_bound(frequency.begin(),
                                                 frequency.end(),
                                                 max(0.0, root - halfWindow));
            const auto windowUpper =
                upper_bound(frequency.begin(),
                            frequency.end(),
                            min(maxFrequency, root + halfWindow));
            const size_t windowBegin =
                max<size_t>(0, distance(frequency.begin(), windowLower));
            const size_t windowEnd = min<size_t>(
                frequency.size() - 1, distance(frequency.begin(), windowUpper));
            const bool useNarrowLimit =
                width < options.narrowPeakSpacingFactor * spacing
                || width < options.narrowPeakOmegaFactor * max(1.0, root);
            for (size_t k = 0; k < kernels.size(); ++k) {
              const auto &kernel = kernels[k];
              const double baseline =
                  integrateBaseline(kernel,
                                    frequency[windowBegin],
                                    frequency[windowEnd]);
              double refined = 1.5 * peakArea * kernel(root);
              if (!useNarrowLimit) {
                const double lower = frequency[windowBegin];
                const double upper = frequency[windowEnd];
                const double uLower = atan((lower - root) / width);
                const double uUpper = atan((upper - root) / width);
                auto transformedPeak = [&](const double &u) {
                  const double cosine = cos(u);
                  const double Omega =
                      clamp(root + width * tan(u), lower, upper);
                  return 1.5
                         * computeDSFValue(
                             x, Omega, Theta, mu, interaction, responseItg)
                         * kernel(Omega) * width / (cosine * cosine);
                };
                peakItg.compute(transformedPeak, ItgParam(uLower, uUpper));
                refined = peakItg.getSolution();
              }
              const double errorScale =
                  max({abs(refined), abs(baseline), 1.0e-10});
              if (abs(refined - baseline) > in->getIntError() * errorScale) {
                result(i, k) += refined - baseline;
              }
            }
          }
        }
        if (right == scanLimit) break;
        left = right;
        leftValue = rightValue;
      }
    }
    return result;
  }

} // namespace thermoUtil
