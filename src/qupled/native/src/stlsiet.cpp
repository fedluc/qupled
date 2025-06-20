#include "stlsiet.hpp"
#include "format.hpp"
#include "mpi_util.hpp"

using namespace std;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

StlsIet::StlsIet(const std::shared_ptr<const StlsIetInput> &in_)
    : Stls(in_, true),
      iet(in_, wvg),
      itg2D(std::make_shared<Integrator2D>(ItgType::DEFAULT,
                                           in_->getIntError())) {
  const bool segregatedItg = in().getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = wvg; }
}

void StlsIet::init() {
  Stls::init();
  iet.init();
}

void StlsIet::computeLfc() {
  // Setup interpolators
  const shared_ptr<Interpolator1D> ssfItp =
      make_shared<Interpolator1D>(wvg, ssf);
  const shared_ptr<Interpolator1D> lfcItp =
      make_shared<Interpolator1D>(wvg[0], lfc(0, 0), wvg.size());
  const shared_ptr<Interpolator1D> bfItp =
      make_shared<Interpolator1D>(wvg, getBf());
  // Compute the stls constribution to the local field correction
  Stls::computeLfc();
  // Compute the iet contribution to the local field correction
  for (size_t i = 0; i < wvg.size(); ++i) {
    StlsIetUtil::Slfc lfcTmp(
        wvg[i], wvg.front(), wvg.back(), ssfItp, lfcItp, bfItp, itgGrid, itg2D);
    lfc(i, 0) += lfcTmp.get();
  }
}

bool StlsIet::initialGuessFromInput() {
  const bool ssfIsSetFromInput = Stls::initialGuessFromInput();
  if (!ssfIsSetFromInput) { return false; }
  const bool lfcIsSetFromInput = iet.initialGuessFromInput(lfc);
  if (!lfcIsSetFromInput) { return false; }
  return true;
}

// -----------------------------------------------------------------
// Slfc class
// -----------------------------------------------------------------

// Compute static local field correction from interpolator
double StlsIetUtil::Slfc::lfc(const double &y) const { return lfci->eval(y); }

// Compute bridge function from interpolator
double StlsIetUtil::Slfc::bf(const double &y) const { return bfi->eval(y); }

// Get at finite temperature
double StlsIetUtil::Slfc::get() const {
  if (x == 0.0) return 0.0;
  auto wMin = [&](const double &y) -> double {
    return (y > x) ? y - x : x - y;
  };
  auto wMax = [&](const double &y) -> double { return min(yMax, x + y); };
  auto func1 = [&](const double &y) -> double { return integrand1(y); };
  auto func2 = [&](const double &w) -> double { return integrand2(w); };
  itg->compute(func1, func2, Itg2DParam(yMin, yMax, wMin, wMax), itgGrid);
  return 3.0 / (8.0 * x) * itg->getSolution() + bf(x);
}

// Level 1 integrand
double StlsIetUtil::Slfc::integrand1(const double &y) const {
  if (y == 0.0) return 0.0;
  return (-bf(y) - (ssf(y) - 1.0) * (lfc(y) - 1.0)) / y;
}

// Level 2 integrand
double StlsIetUtil::Slfc::integrand2(const double &w) const {
  const double y = itg->getX();
  const double y2 = y * y;
  const double w2 = w * w;
  const double x2 = x * x;
  return (w2 - y2 - x2) * w * (ssf(w) - 1.0);
}