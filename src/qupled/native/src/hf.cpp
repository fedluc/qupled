#include "hf.hpp"
#include "chemical_potential.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"

using namespace std;
using namespace thermoUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;

// Constructor
HF::HF(const Input &in_, const bool verbose_)
    : Logger(verbose_ && isRoot()),
      in(in_),
      itg(std::make_shared<Integrator1D>(ItgType::DEFAULT, in_.getIntError())) {
  // Allocate arrays to the correct size
  const size_t nx = in.getWaveVectorGrid().size();
  const size_t nl = in.getNMatsubara();
  idr.resize(nx, nl);
  slfc.resize(nx);
  ssf.resize(nx);
}

// Compute scheme
int HF::compute() {
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
void HF::init() {
  print("Computing chemical potential: ");
  computeChemicalPotential();
  println("Done");
  print("Computing ideal density response: ");
  computeIdr();
  println("Done");
}

// Compute chemical potential
void HF::computeChemicalPotential() {
  if (in.getDegeneracy() == 0.0) return;
  const vector<double> &guess = in.getChemicalPotentialGuess();
  ChemicalPotential mu_(in.getDegeneracy());
  mu_.compute(guess);
  mu = mu_.get();
}

// Compute ideal density response
void HF::computeIdr() {
  (in.getDegeneracy() == 0.0) ? computeIdrGround() : computeIdrFinite();
}

void HF::computeIdrFinite() {
  const size_t nx = idr.size(0);
  const size_t nl = idr.size(1);
  const vector<double> &wvg = in.getWaveVectorGrid();
  for (size_t i = 0; i < nx; ++i) {
    HFUtil::Idr idrTmp(
        nl, wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    idr.fill(i, idrTmp.get());
  }
}

void HF::computeIdrGround() {
  const size_t nx = idr.size(0);
  const size_t nl = idr.size(1);
  const vector<double> &wvg = in.getWaveVectorGrid();
  for (size_t i = 0; i < nx; ++i) {
    for (size_t l = 0; l < nl; ++l) {
      HFUtil::IdrGround idrTmp(wvg[i], l);
      idr(i, l) = idrTmp.get();
    }
  }
}

// Compute Static structure factor
void HF::computeSsf() {
  (in.getDegeneracy() == 0.0) ? computeSsfGround() : computeSsfFinite();
}

void HF::computeSsfFinite() {
  const vector<double> &wvg = in.getWaveVectorGrid();
  for (size_t i = 0; i < wvg.size(); ++i) {
    HFUtil::Ssf ssfTmp(
        wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    ssf[i] = ssfTmp.get();
  }
}

void HF::computeSsfGround() {
  const vector<double> &wvg = in.getWaveVectorGrid();
  for (size_t i = 0; i < wvg.size(); ++i) {
    HFUtil::SsfGround ssfTmp(wvg[i]);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void HF::computeSlfc() {
  assert(slfc.size() == in.getWaveVectorGrid().size());
  for (auto &s : slfc) {
    s = 1;
  }
}

// Getters
vector<double> HF::getRdf(const vector<double> &r) const {
  const vector<double> &wvg = in.getWaveVectorGrid();
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the radial distribution function");
    return vector<double>();
  }
  return computeRdf(r, wvg, ssf);
}

vector<double> HF::getSdr() const {
  const double theta = in.getDegeneracy();
  if (isnan(theta) || theta == 0.0) { return vector<double>(); }
  const vector<double> &wvg = in.getWaveVectorGrid();
  vector<double> sdr(wvg.size(), -1.5 * theta);
  const double fact = 4 * lambda * in.getCoupling() / M_PI;
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x2 = wvg[i] * wvg[i];
    const double phi0 = idr(i, 0);
    sdr[i] *= phi0 / (1.0 + fact / x2 * (1.0 - slfc[i]) * phi0);
  }
  return sdr;
}

double HF::getUInt() const {
  const vector<double> &wvg = in.getWaveVectorGrid();
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
double HFUtil::Idr::integrand(const double &y, const int &l) const {
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
double HFUtil::Idr::integrand(const double &y) const {
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
vector<double> HFUtil::Idr::get() const {
  assert(Theta > 0.0);
  vector<double> res(nl);
  const auto itgParam = ItgParam(yMin, yMax);
  for (int l = 0; l < nl; ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand(y) : integrand(y, l);
    };
    itg->compute(func, itgParam);
    res[l] = itg->getSolution();
  }
  return res;
}

// -----------------------------------------------------------------
// IdrGround class
// -----------------------------------------------------------------

// Get
double HFUtil::IdrGround::get() const {
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
double HFUtil::Ssf::integrand(const double &y) const {
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
double HFUtil::Ssf::get() const {
  assert(Theta > 0.0);
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(yMin, yMax));
  return 1.0 + itg->getSolution();
}

// -----------------------------------------------------------------
// SsfHFGround class
// -----------------------------------------------------------------

// Static structure factor at zero temperature
double HFUtil::SsfGround::get() const {
  if (x < 2.0) {
    return (x / 16.0) * (12.0 - x * x);
  } else {
    return 1.0;
  }
}