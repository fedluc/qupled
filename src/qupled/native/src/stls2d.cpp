#include "stls2d.hpp"
#include "stls.hpp"
#include "bin_util.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <fmt/core.h>
#include <gsl/gsl_sf_ellint.h>

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
    computeSsf2D();
    // Update static local field correction
    computeSlfc2D();
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

void Stls2D::computeChemicalPotential2D() {
  mu = log(exp(1/in.getDegeneracy())-1);
}

void Stls2D::computeIdr2D() {
  if (in.getDegeneracy() == 0.0) return;
  const size_t nx = wvg.size();
  const size_t nl = in.getNMatsubara();
  assert(idr.size(0) == nx && idr.size(1) == nl);
  for (size_t i = 0; i < nx; ++i) {
    Idr idrTmp(
        nl, wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    idr.fill(i, idrTmp.get());
  }
}

void Stls2D::computeSsf2D() {
  // placeholder;
}

void Stls2D::computeSsfHF2D() {
  Integrator2D itg2(ItgType::DEFAULT, ItgType::DEFAULT, in.getIntError());
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  for (size_t i = 0; i < wvg.size(); ++i) {
    SSFHF2D ssfHF2DTmp(in.getDegeneracy(),
                        wvg[i],
                        mu,
                        wvg.front(),
                        wvg.back(),
                        itgGrid,
                        itg2);
    ssfHF[i] = ssfHF2DTmp.get();
  }
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

// -----------------------------------------------------------------
// Idr 2D class
// -----------------------------------------------------------------

vector<double> Idr::get2DStls() const {
  assert(Theta > 0.0);
  vector<double> res(nl);
  for (int l = 0; l < nl; ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand2DStls(y) : integrand2DStls(y, l);
    };
    double upperLimit = (l == 0) ? yMax / 2.0 : yMax;
    const auto itgParam = ItgParam(yMin, upperLimit);
    itg.compute(func, itgParam);
    if (l == 0) {
      res[l] = 1.0 - exp(-1.0 / Theta) - itg.getSolution();
    } else {
      res[l] = itg.getSolution();
    }
  }
  return res;
}

// Integrand for frequency = l and wave-vector = x
double Idr::integrand2DStls(const double &y, const int &l) const {
  double tphi;
  double y2 = y * y;
  double x2 = x * x;
  double x4 = x2 * x2;
  double plT = M_PI * l * Theta;
  double plT2 = plT * plT;
  double exp1 = x4 / 4.0 - x2 * y2 - plT2;
  if (exp1 < 0.0) {
    tphi = atan(x2 * plT / exp1);
  } else {
    tphi = M_PI - atan(x2 * plT / exp1);
  }
  if (x > 0.0) {
    return y / (exp(y2 / Theta - mu) + 1.0)
           * 2.0 * abs(cos(tphi))/ pow((exp1 * exp1 + x4 * plT2), 0.25);
  } else {
    return 0;
  }
}

// Integrand for frequency = 0 and vector = x
double Idr::integrand2DStls(const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  if (x > 0.0) {
    - 1.0 / (Theta * x * pow(cosh(y2 / (2 * Theta) - mu/2), 2) * y * sqrt(x2 / 4.0 - y2));
  } else {
    return 0; 
  }
}

// -----------------------------------------------------------------
// SLFC 2D STLS class
// -----------------------------------------------------------------

// Integrand 2D STLS
double Slfc::integrand2DStls(const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  double xmy = (x - y) / (x * M_PI);
  double xpy = (x + y) / (x * M_PI);
  double argElli = 2 * sqrt(x * y) / (x + y);
  if (x > 0.0) {
  return - (ssf(y) - 1.0) * (xmy * Integrator1D::ellipticK(argElli) + xpy * Integrator1D::ellipticE(argElli));}
  else {
    return 0;
  }
}

// Get result of integration
double Slfc::get2DStls() const {
  auto func = [&](const double &y) -> double { return integrand2DStls(y); };
  itg.compute(func, ItgParam(yMin, yMax));
  return itg.getSolution();
}

// -----------------------------------------------------------------
// Elliptic function integrators
// -----------------------------------------------------------------

double Integrator1D::ellipticK(const double &k) {
  gsl_sf_result result;
  int status = gsl_sf_ellint_Kcomp_e(k, GSL_PREC_DOUBLE, &result);
  if (status != GSL_SUCCESS) {
    throw std::runtime_error("GSL error in ellipticK");
  }
  return result.val;
}

double Integrator1D::ellipticE(const double &k) {
  gsl_sf_result result;
  int status = gsl_sf_ellint_Ecomp_e(k, GSL_PREC_DOUBLE, &result);
  if (status != GSL_SUCCESS) {
    throw std::runtime_error("GSL error in ellipticE");
  }
  return result.val;
}

// -----------------------------------------------------------------
// SSF HF 2D
// -----------------------------------------------------------------

inline double coth(double x) {
  return 1.0 / tanh(x);
}

// Outer integrand
double SSFHF2D::integrandOut(const double y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / Theta - mu) * M_PI + M_PI);
}

// Inner integrand
double SSFHF2D::integrandIn(const double p) const {
  const double y = itg2.getX();
  double x2 = x * x;
  return coth(x2 / (2 * Theta) - x * y / Theta * cos(p));
}

// Get total QAdder
double SSFHF2D::get() const {
  auto func1 = [&](const double &y) -> double {
    return integrandOut(y);
  };
  auto func2 = [&](const double &p) -> double {
    return integrandIn(p);
  };
  itg2.compute(
      func1,
      func2,
      Itg2DParam(limits.first, limits.second, 0, 2 * M_PI),
      itgGrid);
  return itg2.getSolution();
}