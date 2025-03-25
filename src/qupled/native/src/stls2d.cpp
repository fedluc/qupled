#include "stls2d.hpp"
#include "stls.hpp"
#include "bin_util.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <fmt/core.h>

using namespace std;
using namespace vecUtil;
using namespace binUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

Stls2D::Stls2D(const StlsInput &in_, const bool verbose_, const bool writeFiles_)
    : Stls(in_, verbose_, writeFiles_),
      in(in_),
      writeFiles(writeFiles_ && isRoot()) {
  // Set name of recovery files
  recoveryFileName = fmt::format("recovery_rs{:.3f}_theta{:.3f}_{}.bin",
                                 in.getCoupling(),
                                 in.getDegeneracy(),
                                 in.getTheory());
  // Allocate arrays
  const size_t nx = wvg.size();
  slfcNew.resize(nx);
  if (useIet) { bf.resize(nx); }
}

int Stls2D::compute() {
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

// Initialize basic properties
void Stls2D::init() {
  print("Computing 2D chemical potential: ");
  computeChemicalPotential2D();
  println("Done");
  print("Computing 2D ideal density response: ");
  computeIdr2D();
  println("Done");
  print("Computing 2D Hartree-Fock static structure factor: ");
  computeSsfHF2D();
  println("Done");
}
void Stls2D::computeChemicalPotential2D() {
  // placeholder
}

void Stls2D::computeIdr2D() {
  // placeholder
}

void Stls2D::computeSsfHF2D() {
  // placeholder
}

// Compute static local field correction
void Stls2D::computeSlfc2D() {
  assert(ssf.size() == wvg.size());
  assert(slfc.size() == wvg.size());
  computeSlfcStls2D();
}

void Stls2D::computeSlfcStls2D() {
  const int nx = wvg.size();
  const Interpolator1D itp(wvg, ssf);
  for (int i = 0; i < nx; ++i) {
    Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itp, itg);
    slfcNew[i] = slfcTmp.get();
  }
}

// stls iterations
void Stls2D::doIterations() {
  const int maxIter = in.getNIter();
  const int outIter = in.getOutIter();
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
    // Write output
    if (counter % outIter == 0 && writeFiles) { writeRecovery(); }
    // End timing
    double toc = timer();
    // Print diagnostic
    println(fmt::format("--- iteration {:d} ---", counter));
    println(fmt::format("Elapsed time: {:.3f} seconds", toc - tic));
    println(fmt::format("Residual error: {:.5e}", err));
    fflush(stdout);
  }
}


// -----------------------------------------------------------------
// SlfcBase class
// -----------------------------------------------------------------

// Compute static structure factor from interpolator
double SlfcBase::ssf(const double &y) const { return ssfi.eval(y); }

// -----------------------------------------------------------------
// Slfc class
// -----------------------------------------------------------------

// Get result of integration
double Slfc::get() const {
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(yMin, yMax));
  return itg.getSolution();
}

// Integrand
double Slfc::integrand(const double &y) const {
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
