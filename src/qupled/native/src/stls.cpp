#include "stls.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <format>
#include <sstream>

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

Stls::Stls(const StlsInput &in_, const bool verbose_)
    : Rpa(in_, verbose_),
      in(in_) {
  // Allocate arrays
  const size_t nx = wvg.size();
  slfcNew.resize(nx);
}

int Stls::compute() {
  try {
    init();
    println("Structural properties calculation ...");
    doIterations();
    println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

// Compute static local field correction
void Stls::computeSlfc() {
  assert(ssf.size() == wvg.size());
  assert(slfc.size() == wvg.size());
  const int nx = wvg.size();
  const shared_ptr<Interpolator1D> itp = make_shared<Interpolator1D>(wvg, ssf);
  for (int i = 0; i < nx; ++i) {
    StlsUtil::Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itp, itg);
    slfcNew[i] = slfcTmp.get();
  }
}

// stls iterations
void Stls::doIterations() {
  const int maxIter = in.getNIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  while (counter < maxIter + 1 && err > minErr) {
    // Start timing
    double tic = timer();
    // Update static structure factor
    computeSsf();
    // Update static local field correction
    computeSlfc();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // End timing
    double toc = timer();
    // Print diagnostic
    println(format("--- iteration {:d} ---", counter));
    println(format("Elapsed time: {:.3f} seconds", toc - tic));
    println(format("Residual error: {:.5e}", err));
    fflush(stdout);
  }
}

// Initial guess for stls iterations
void Stls::initialGuess() {
  // From guess in input
  if (initialGuessFromInput()) { return; }
  // Default
  fill(slfc.begin(), slfc.end(), 0.0);
}

bool Stls::initialGuessFromInput() {
  const Interpolator1D slfci(in.getGuess().wvg, in.getGuess().slfc);
  if (!slfci.isValid()) { return false; }
  const double xmaxi = in.getGuess().wvg.back();
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x = wvg[i];
    if (x <= xmaxi) {
      slfc[i] = slfci.eval(x);
    } else {
      slfc[i] = 1.0;
    }
  }
  return true;
}

// Compute residual error for the stls iterations
double Stls::computeError() const { return rms(slfc, slfcNew, false); }

// Update solution during stls iterations
void Stls::updateSolution() {
  const double aMix = in.getMixingParameter();
  slfc = linearCombination(slfcNew, aMix, slfc, 1 - aMix);
}

// -----------------------------------------------------------------
// SlfcBase class
// -----------------------------------------------------------------

// Compute static structure factor from interpolator
double StlsUtil::SlfcBase::ssf(const double &y) const { return ssfi->eval(y); }

// -----------------------------------------------------------------
// Slfc class
// -----------------------------------------------------------------

// Get result of integration
double StlsUtil::Slfc::get() const {
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(yMin, yMax));
  return itg->getSolution();
}

// Integrand
double StlsUtil::Slfc::integrand(const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  if (x == 0.0 || y == 0.0) { return 0.0; }
  if (x == y) { return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0); };
  if (x > y) {
    return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0)
           * (1 + (x2 - y2) / (2 * x * y) * log((x + y) / (x - y)));
  }
  return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0)
         * (1 + (x2 - y2) / (2 * x * y) * log((x + y) / (y - x)));
}