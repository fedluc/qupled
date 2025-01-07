#include "rpa.hpp"
#include "cdual.hpp"
#include "chemical_potential.hpp"
#include "dual.hpp"
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
  const size_t nx = wvg.size();
  double wpGuess = 4.0 * sqrt(lambda * rs / 3.0 / M_PI); // Plasma frequency
  assert(slfc.size() == nx);
  assert(ssf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    const double x = wvg[i];
    DielectricResponse dr(x, rs, slfc[i]);
    const double wp = (wpGuess >= 0) ? dr.plasmon(wpGuess) : -1;
    wpGuess = wp;
    SsfGround ssfTmp(x, rs, slfc[i], wp, itg);
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
  if (isnan(theta)) { return vector<double>(); }
  if (theta == 0.0) {
    println("The static density response cannot be computed in the "
            "ground state.");
    return vector<double>();
  }
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

// Real part at zero temperature
template <typename T>
T IdrGround::real() const {
  const T dOmega = T(Omega, 0);
  T adder1 = T(0.0);
  T adder2 = T(0.0);
  T preFactor = T(0.0);
  if (x > 0.0) {
    const double x_2 = x / 2.0;
    const T Omega_2x = dOmega / (2.0 * x);
    const T sumFactor = x_2 + Omega_2x;
    const T diffFactor = x_2 - Omega_2x;
    const T sumFactor2 = sumFactor * sumFactor;
    const T diffFactor2 = diffFactor * diffFactor;
    preFactor = preFactor + 0.5;
    if (sumFactor.val() != 1.0) {
      T log_sum_arg = (sumFactor + 1.0) / (sumFactor - 1.0);
      if (log_sum_arg.val() < 0.0) log_sum_arg = -1.0 * log_sum_arg;
      adder1 = 1.0 / (4.0 * x) * (1.0 - sumFactor2) * log(log_sum_arg);
    }
    if (diffFactor.val() != 1.0 && diffFactor.val() != -1.0) {
      T log_diff_arg = (diffFactor + 1.0) / (diffFactor - 1.0);
      if (log_diff_arg.val() < 0.0) log_diff_arg = -1.0 * log_diff_arg;
      adder2 = 1.0 / (4.0 * x) * (1.0 - diffFactor2) * log(log_diff_arg);
    }
  }
  return preFactor + adder1 + adder2;
}
template Dual0 IdrGround::real<Dual0>() const;
template Dual11 IdrGround::real<Dual11>() const;
template Dual21 IdrGround::real<Dual21>() const;

// Imaginary part at zero temperature
template <typename T>
T IdrGround::imag() const {
  const T dOmega = T(Omega, 0);
  T preFactor = T(0.0);
  T adder1 = T(0.0);
  T adder2 = T(0.0);
  if (x > 0.0) {
    T x_2 = T(x / 2.0);
    T Omega_2x = dOmega / (2.0 * x);
    T sumFactor = x_2 + Omega_2x;
    T diffFactor = x_2 - Omega_2x;
    T sumFactor2 = sumFactor * sumFactor;
    T diffFactor2 = diffFactor * diffFactor;
    preFactor = T(-M_PI / (4.0 * x));
    if (sumFactor2.val() < 1.0) { adder1 = 1 - sumFactor2; }
    if (diffFactor2.val() < 1.0) { adder2 = 1 - diffFactor2; }
  }
  return preFactor * (adder1 - adder2);
}
template Dual0 IdrGround::imag<Dual0>() const;
template Dual11 IdrGround::imag<Dual11>() const;
template Dual21 IdrGround::imag<Dual21>() const;

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
  const double fact1 = 4.0 * lambda * rs / M_PI;
  const double x2 = x * x;
  double fact2 = 0.0;
  for (int l = 0; l < nl; ++l) {
    const double fact3 = 1.0 + fact1 / x2 * (1 - slfc) * idr[l];
    double fact4 = idr[l] * idr[l] / fact3;
    if (l > 0) fact4 *= 2;
    fact2 += fact4;
  }
  return ssfHF - 1.5 * fact1 / x2 * Theta * (1 - slfc) * fact2;
}

// -----------------------------------------------------------------
// SsfGround class
// -----------------------------------------------------------------

// Get result of integration
double SsfGround::get() const {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(yMin, yMax));
  const double ssfP = plasmon();
  return 1.5 / (M_PI)*itg.getSolution() + ssfP;
}

// Integrand for zero temperature calculations
double SsfGround::integrand(const double &Omega) const {
  const DielectricResponse dr = DielectricResponse(x, rs, slfc);
  const CDual0 phi = 1.0 / dr.get<Dual0>(Omega);
  return -phi.imag.val();
}

// Plasmon contribution to the static structure factor
double SsfGround::plasmon() const {
  if (wp < 0) { return 0.0; }
  const DielectricResponse dr = DielectricResponse(x, rs, slfc);
  return 1.5 / abs(dr.get<Dual11>(wp).real.dx());
}

// -----------------------------------------------------------------
// DielectricResponse class
// -----------------------------------------------------------------

// Real part and its derivative
template <typename T>
CDual<T> DielectricResponse::get(const double &Omega) const {
  const IdrGround idr = IdrGround(Omega, x);
  const CDual<T> cidr = CDual<T>(idr.real<T>(), idr.imag<T>());
  return 1.0 + cidr / (1.0 + cidr * (ip * (1 - slfc) - 1.0));
}

// Get the plasmon frequency
double DielectricResponse::plasmon(const double &guess) const {
  if (x == 0.0) { return wp; }
  // Compute plasmon frequency
  auto func = [this](const double &Omega) -> pair<double, double> {
    const Dual21 deq = dispersionEquation(Omega);
    return pair<double, double>(deq.dx(), deq.dxx());
  };
  QuasiNewtonRootSolver rsol;
  try {
    rsol.solve(func, guess);
  } catch (const std::runtime_error &e) {
    // The plasmon does not exist
    return -1;
  }
  // Check if the minimum is a zero
  const double wp = rsol.getSolution();
  bool isZero = abs(dispersionEquation(wp).val()) < 1e-10;
  // Output
  return (isZero) ? wp : -1;
}

// Dispersion equation
Dual21 DielectricResponse::dispersionEquation(const double &Omega) const {
  const IdrGround idr = IdrGround(Omega, x);
  const CDual21 cidr = CDual21(idr.real<Dual21>(), idr.imag<Dual21>());
  const CDual21 deq = 1.0 + ip * cidr * (1.0 - slfc);
  return deq.real * deq.real + deq.imag * deq.imag;
}
