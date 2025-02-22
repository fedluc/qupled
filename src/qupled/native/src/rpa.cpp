#include "rpa.hpp"
#include "chemical_potential.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include <cmath>

using namespace std;
using namespace thermoUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;

// Constructor
Rpa::Rpa(const RpaInput &in_, const bool verbose_)
    : Logger(verbose_ && isRoot()),
      in(in_),
      itg(ItgType::DEFAULT, in_.getIntError()) {
  // Assemble the wave-vector grid
  buildWvGrid();
  // Allocate arrays to the correct size
  const size_t nx = wvg.size();
  const size_t nl = in.getNMatsubara();
  idr.resize(nx, nl);
  slfc.resize(nx);
  ssf.resize(nx);
  ssfHF.resize(nx);
}

// Compute scheme
int Rpa::compute() {
  try {
    init();
    println("Structural properties calculation ...");
    print("Computing static local field correction: ");
    computeSlfc();
    println("Done");
    print("Computing static structure factor: ");
    computeSsf();
    println("Done");
    println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

// Initialize basic properties
void Rpa::init() {
  print("Computing chemical potential: ");
  computeChemicalPotential();
  println("Done");
  print("Computing ideal density response: ");
  computeIdr();
  println("Done");
  print("Computing Hartree-Fock static structure factor: ");
  computeSsfHF();
  println("Done");
}

// Set up wave-vector grid
void Rpa::buildWvGrid() {
  wvg.push_back(0.0);
  const double dx = in.getWaveVectorGridRes();
  const double xmax = in.getWaveVectorGridCutoff();
  if (xmax < dx) {
    throwError(
        "The wave-vector grid cutoff must be larger than the resolution");
  }
  while (wvg.back() < xmax) {
    wvg.push_back(wvg.back() + dx);
  }
}

// Compute chemical potential
void Rpa::computeChemicalPotential() {
  if (in.getDegeneracy() == 0.0) return;
  const vector<double> &guess = in.getChemicalPotentialGuess();
  ChemicalPotential mu_(in.getDegeneracy());
  mu_.compute(guess);
  mu = mu_.get();
}

// Compute ideal density response
void Rpa::computeIdr() {
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

// Compute Hartree-Fock static structure factor
void Rpa::computeSsfHF() {
  assert(ssfHF.size() == wvg.size());
  if (in.getDegeneracy() == 0.0) {
    computeSsfHFGround();
  } else {
    computeSsfHFFinite();
  }
}

void Rpa::computeSsfHFFinite() {
  for (size_t i = 0; i < wvg.size(); ++i) {
    SsfHF ssfTmp(wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    ssfHF[i] = ssfTmp.get();
  }
}

void Rpa::computeSsfHFGround() {
  for (size_t i = 0; i < wvg.size(); ++i) {
    SsfHFGround ssfTmp(wvg[i]);
    ssfHF[i] = ssfTmp.get();
  }
}

// Compute static structure factor
void Rpa::computeSsf() {
  assert(ssf.size() == wvg.size());
  if (in.getDegeneracy() == 0.0) {
    computeSsfGround();
  } else {
    computeSsfFinite();
  }
}

// Compute static structure factor at finite temperature
void Rpa::computeSsfFinite() {
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const size_t nx = wvg.size();
  const size_t nl = idr.size(1);
  assert(slfc.size() == nx);
  assert(ssf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], slfc[i], nl, &idr(i));
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Rpa::computeSsfGround() {
  const double rs = in.getCoupling();
  const double OmegaMax = in.getFrequencyCutoff();
  const size_t nx = wvg.size();
  assert(slfc.size() == nx);
  assert(ssf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    const double x = wvg[i];
    SsfGround ssfTmp(x, rs, ssfHF[i], slfc[i], OmegaMax, itg);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void Rpa::computeSlfc() {
  assert(slfc.size() == wvg.size());
  for (auto &s : slfc) {
    s = 0;
  }
}

// Getters
vector<double> Rpa::getRdf(const vector<double> &r) const {
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the radial distribution function");
    return vector<double>();
  }
  return computeRdf(r, wvg, ssf);
}

vector<double> Rpa::getSdr() const {
  const double theta = in.getDegeneracy();
  if (isnan(theta) || theta == 0.0) { return vector<double>(); }
  vector<double> sdr(wvg.size(), -1.5 * theta);
  const double fact = 4 * lambda * in.getCoupling() / M_PI;
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x2 = wvg[i] * wvg[i];
    const double phi0 = idr(i, 0);
    sdr[i] *= phi0 / (1.0 + fact / x2 * (1.0 - slfc[i]) * phi0);
  }
  return sdr;
}

double Rpa::getUInt() const {
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the internal energy");
    return numUtil::Inf;
  }
  return computeInternalEnergy(wvg, ssf, in.getCoupling());
}

// -----------------------------------------------------------------
// Idr class
// -----------------------------------------------------------------

// Integrand for frequency = l and wave-vector = x
double Idr::integrand(const double &y, const int &l) const {
  double y2 = y * y;
  double x2 = x * x;
  double txy = 2 * x * y;
  double tplT = 2 * M_PI * l * Theta;
  double tplT2 = tplT * tplT;
  if (x > 0.0) {
    return 1.0 / (2 * x) * y / (exp(y2 / Theta - mu) + 1.0)
           * log(((x2 + txy) * (x2 + txy) + tplT2)
                 / ((x2 - txy) * (x2 - txy) + tplT2));
  } else {
    return 0;
  }
}

// Integrand for frequency = 0 and vector = x
double Idr::integrand(const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  double xy = x * y;
  if (x > 0.0) {
    if (x < 2 * y) {
      return 1.0 / (Theta * x)
             * ((y2 - x2 / 4.0) * log((2 * y + x) / (2 * y - x)) + xy) * y
             / (exp(y2 / Theta - mu) + exp(-y2 / Theta + mu) + 2.0);
    } else if (x > 2 * y) {
      return 1.0 / (Theta * x)
             * ((y2 - x2 / 4.0) * log((2 * y + x) / (x - 2 * y)) + xy) * y
             / (exp(y2 / Theta - mu) + exp(-y2 / Theta + mu) + 2.0);
    } else {
      return 1.0 / (Theta)*y2
             / (exp(y2 / Theta - mu) + exp(-y2 / Theta + mu) + 2.0);
      ;
    }
  } else {
    return (2.0 / Theta) * y2
           / (exp(y2 / Theta - mu) + exp(-y2 / Theta + mu) + 2.0);
  }
}

// Get result of integration
vector<double> Idr::get() const {
  assert(Theta > 0.0);
  vector<double> res(nl);
  const auto itgParam = ItgParam(yMin, yMax);
  for (int l = 0; l < nl; ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand(y) : integrand(y, l);
    };
    itg.compute(func, itgParam);
    res[l] = itg.getSolution();
  }
  return res;
}

// -----------------------------------------------------------------
// IdrGround class
// -----------------------------------------------------------------

// Get
double IdrGround::get() const {
  const double x2 = x * x;
  const double Omega2 = Omega * Omega;
  const double tx = 2.0 * x;
  const double x2ptx = x2 + tx;
  const double x2mtx = x2 - tx;
  const double x2ptx2 = x2ptx * x2ptx;
  const double x2mtx2 = x2mtx * x2mtx;
  const double logarg = (x2ptx2 + Omega2) / (x2mtx2 + Omega2);
  const double part1 = (0.5 - x2 / 8.0 + Omega2 / (8.0 * x2)) * log(logarg);
  const double part2 =
      0.5 * Omega * (atan(x2ptx / Omega) - atan(x2mtx / Omega));
  if (x > 0.0) { return (part1 - part2 + x) / tx; }
  return 0;
}

// -----------------------------------------------------------------
// SsfHF class
// -----------------------------------------------------------------

// Integrand
double SsfHF::integrand(const double &y) const {
  double y2 = y * y;
  double ypx = y + x;
  double ymx = y - x;
  if (x > 0.0) {
    return -3.0 * Theta / (4.0 * x) * y / (exp(y2 / Theta - mu) + 1.0)
           * log((1 + exp(mu - ymx * ymx / Theta))
                 / (1 + exp(mu - ypx * ypx / Theta)));
  } else {
    return -3.0 * y2
           / ((1.0 + exp(y2 / Theta - mu)) * (1.0 + exp(y2 / Theta - mu)));
  }
}

// Get result of integration
double SsfHF::get() const {
  assert(Theta > 0.0);
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(yMin, yMax));
  return 1.0 + itg.getSolution();
}

// -----------------------------------------------------------------
// SsfHFGround class
// -----------------------------------------------------------------

// Static structure factor at zero temperature
double SsfHFGround::get() const {
  if (x < 2.0) {
    return (x / 16.0) * (12.0 - x * x);
  } else {
    return 1.0;
  }
}

// -----------------------------------------------------------------
// Ssf class
// -----------------------------------------------------------------

// Get at finite temperature for any scheme
double Ssf::get() const {
  assert(Theta > 0.0);
  if (rs == 0.0) return ssfHF;
  if (x == 0.0) return 0.0;
  double fact2 = 0.0;
  for (int l = 0; l < nl; ++l) {
    const double fact3 = 1.0 + ip * (1 - slfc) * idr[l];
    double fact4 = idr[l] * idr[l] / fact3;
    if (l > 0) fact4 *= 2;
    fact2 += fact4;
  }
  return ssfHF - 1.5 * ip * Theta * (1 - slfc) * fact2;
}

// -----------------------------------------------------------------
// SsfGround class
// -----------------------------------------------------------------

double SsfGround::get() {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(0, OmegaMax));
  return 1.5 / (M_PI)*itg.getSolution() + ssfHF;
}

double SsfGround::integrand(const double &Omega) const {
  const double idr = IdrGround(x, Omega).get();
  return idr / (1.0 + ip * idr * (1.0 - slfc)) - idr;
}
