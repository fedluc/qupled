#include "thermo/dsf_moments.hpp"
#include "schemes/input.hpp"
#include "thermo/dsf_integration.hpp"
#include "util/mpi_util.hpp"
#include <cmath>
#include <functional>

using namespace std;

namespace thermoUtil {

  vector<double> computeAdaptiveFirstMomentFromDSF(
      const shared_ptr<const Input> &in,
      const vector<double> &wvg,
      const vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf) {
    const double Theta = in->getDegeneracy();
    if (Theta <= 0.0) {
      MPIUtil::throwError(
          "Adaptive DSF integration requires finite temperature");
    }
    const double beta = 1.0 / Theta;
    vector<function<double(double)>> kernels = {[=](const double &Omega) {
      return Omega * -expm1(-beta * Omega);
    }};
    const Vector2D result = computeAdaptiveWeightedDSFIntegral(
        in, wvg, frequency, mu, lfc, dsf, kernels, 2.0);
    vector<double> moment(wvg.size());
    for (size_t i = 0; i < wvg.size(); ++i)
      moment[i] = result(i, 0);
    return moment;
  }

  vector<double> computeAdaptiveThirdMomentFromDSF(
      const shared_ptr<const Input> &in,
      const vector<double> &wvg,
      const vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf) {
    const double Theta = in->getDegeneracy();
    if (Theta <= 0.0) {
      MPIUtil::throwError(
          "Adaptive DSF integration requires finite temperature");
    }
    const double beta = 1.0 / Theta;
    vector<function<double(double)>> kernels = {[=](const double &Omega) {
      return Omega * Omega * Omega * -expm1(-beta * Omega);
    }};
    const Vector2D result = computeAdaptiveWeightedDSFIntegral(
        in, wvg, frequency, mu, lfc, dsf, kernels, 1.0);
    vector<double> moment(wvg.size());
    for (size_t i = 0; i < wvg.size(); ++i)
      moment[i] = result(i, 0);
    return moment;
  }

} // namespace thermoUtil
