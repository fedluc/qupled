#include "hf.hpp"
#include "chemical_potential.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"

using namespace std;
using namespace dimensionsUtil;
using namespace thermoUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;
using Itg2DParam = Integrator2D::Param;

// Constructor
HF::HF(const std::shared_ptr<const Input> &in_, const bool verbose_)
    : Logger(verbose_),
      inPtr(std::move(in_)),
      itg(std::make_shared<Integrator1D>(ItgType::DEFAULT,
                                         in_->getIntError())) {
  // Assemble the wave-vector grid
  buildWaveVectorGrid();
  // Allocate arrays to the correct size
  const size_t nx = wvg.size();
  const size_t nl = in().getNMatsubara();
  idr.resize(nx, nl);
  lfc.resize(nx, 1);
  ssf.resize(nx);
}

// Compute scheme
int HF::compute() {
  try {
    init();
    println("Structural properties calculation ...");
    computeStructuralProperties();
    println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void HF::computeStructuralProperties() {
  print("Computing static local field correction: ");
  computeLfc();
  println("Done");

  print("Computing static structure factor: ");
  computeSsf();
  println("Done");
}

void HF::init() {
  if (in().getDegeneracy() == 0.0 && in().getDimension() == Dimension::D2) {
    throw runtime_error("2D Ground state computations are unsupported.");
  }
  print("Computing chemical potential: ");
  computeChemicalPotential();
  println("Done");
  print("Computing ideal density response: ");
  computeIdr();
  println("Done");
}

// Set up wave-vector grid
void HF::buildWaveVectorGrid() {
  wvg.push_back(0.0);
  const double dx = in().getWaveVectorGridRes();
  const double xmax = in().getWaveVectorGridCutoff();
  if (xmax < dx) {
    throwError(
        "The wave-vector grid cutoff must be larger than the resolution");
  }
  while (wvg.back() < xmax) {
    wvg.push_back(wvg.back() + dx);
  }
}

// Compute chemical potential
void HF::computeChemicalPotential() {
  if (in().getDegeneracy() == 0.0) return;
  ChemicalPotential mu_(inPtr);
  mu_.compute(in().getDimension());
  mu = mu_.get();
}

// Compute ideal density response
void HF::computeIdr() {
  (in().getDegeneracy() == 0.0) ? computeIdrGround() : computeIdrFinite();
}

void HF::computeIdrFinite() {
  const size_t nx = idr.size(0);
  for (size_t i = 0; i < nx; ++i) {
    HFUtil::Idr idrTmp(inPtr, wvg[i], mu, wvg.front(), wvg.back(), itg);
    idr.fill(i, idrTmp.get());
  }
}

void HF::computeIdrGround() {
  const size_t nx = idr.size(0);
  const size_t nl = idr.size(1);
  for (size_t i = 0; i < nx; ++i) {
    for (size_t l = 0; l < nl; ++l) {
      HFUtil::IdrGround idrTmp(wvg[i], l);
      idr(i, l) = idrTmp.get();
    }
  }
}

// Compute static structure factor
void HF::computeSsf() {
  (in().getDegeneracy() == 0.0) ? computeSsfGround() : computeSsfFinite();
}

void HF::computeSsfFinite() {
  const bool segregatedItg = in().getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>(in().getIntError());
  for (size_t i = 0; i < wvg.size(); ++i) {
    HFUtil::Ssf ssfTmp(
        inPtr, wvg[i], mu, wvg.front(), wvg.back(), itg, itgGrid, itg2, idr, i);
    ssf[i] = ssfTmp.get();
  }
}

void HF::computeSsfGround() {
  for (size_t i = 0; i < wvg.size(); ++i) {
    HFUtil::SsfGround ssfTmp(wvg[i]);
    ssf[i] = ssfTmp.get();
  }
}

void HF::computeLfc() {
  assert(lfc.size() == wvg.size());
  for (auto &s : lfc) {
    s = 1;
  }
}

// Getters
vector<double> HF::getSdr() const {
  const double theta = in().getDegeneracy();
  const dimensionsUtil::Dimension dim = in().getDimension();
  if (isnan(theta) || theta == 0.0) { return vector<double>(); }
  // Calculate SDR for 2D
  if (dim == Dimension::D2) {
    vector<double> sdr(wvg.size(), -theta);
    const double fact = sqrt(2.0) * in().getCoupling();
    for (size_t i = 0; i < wvg.size(); ++i) {
      const double phi0 = idr(i, 0);
      sdr[i] *= phi0 / (1.0 + fact / wvg[i] * (1.0 - lfc(i, 0)) * phi0);
    }
    return sdr;
  }
  // Calculate SDR for 3D
  vector<double> sdr(wvg.size(), -1.5 * theta);
  const double fact = 4 * numUtil::lambda * in().getCoupling() / M_PI;
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x2 = wvg[i] * wvg[i];
    const double phi0 = idr(i, 0);
    sdr[i] *= phi0 / (1.0 + fact / x2 * (1.0 - lfc(i, 0)) * phi0);
  }
  return sdr;
}

double HF::getUInt() const {
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the internal energy");
    return numUtil::Inf;
  }
  return computeInternalEnergy(
      wvg, ssf, in().getCoupling(), in().getDimension());
}

// -----------------------------------------------------------------
// Idr class
// -----------------------------------------------------------------

std::vector<double> HFUtil::Idr::get() {
  assert(in->getDegeneracy() > 0.0);
  compute(in->getDimension());
  return res;
};

// Compute for 3D systems
void HFUtil::Idr::compute3D() {
  const auto itgParam = ItgParam(yMin, yMax);
  for (int l = 0; l < in->getNMatsubara(); ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand(y) : integrand(y, l);
    };
    itg->compute(func, itgParam);
    res[l] = itg->getSolution();
  }
}

// Compute for 2D systems
void HFUtil::Idr::compute2D() {
  const double Theta = in->getDegeneracy();
  for (int l = 0; l < in->getNMatsubara(); ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand2D(y) : integrand2D(y, l);
    };
    double upperLimit = (l == 0) ? x / 2.0 : yMax;
    const auto itgParam = ItgParam(yMin, upperLimit);
    itg->compute(func, itgParam);
    if (l == 0) {
      res[l] = 1.0 - exp(-1.0 / Theta) - itg->getSolution();
    } else {
      res[l] = itg->getSolution();
    }
  }
}

// Integrand for frequency = l and wave-vector = x
double HFUtil::Idr::integrand(const double &y, const int &l) const {
  const double Theta = in->getDegeneracy();
  const double y2 = y * y;
  const double x2 = x * x;
  const double txy = 2 * x * y;
  const double tplT = 2 * M_PI * l * Theta;
  const double tplT2 = tplT * tplT;
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
  const double Theta = in->getDegeneracy();
  const double y2 = y * y;
  const double x2 = x * x;
  const double xy = x * y;
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

// Integrand for frequency = l and wave-vector = x
double HFUtil::Idr::integrand2D(const double &y, const int &l) const {
  const double Theta = in->getDegeneracy();
  const double y2 = y * y;
  const double x2 = x * x;
  const double x4 = x2 * x2;
  const double plT = M_PI * l * Theta;
  const double plT2 = plT * plT;
  const double exp1 = x4 / 4.0 - x2 * y2 - plT2;
  double phi;
  if (exp1 > 0.0) {
    phi = atan(x2 * plT / exp1) / 2.0;
  } else {
    phi = M_PI / 2.0 - atan(x2 * plT / exp1) / 2.0;
  }
  if (x > 0.0) {
    return y / (exp(y2 / Theta - mu) + 1.0) * 2.0 * abs(cos(phi))
           / pow((exp1 * exp1 + x4 * plT2), 0.25);
  } else {
    return 0;
  }
}

// Integrand for frequency = 0 and vector = x
double HFUtil::Idr::integrand2D(const double &y) const {
  const double Theta = in->getDegeneracy();
  const double y2 = y * y;
  const double x2 = x * x;
  if (x > 0.0) {
    return 1.0 / (Theta * x * pow(cosh(y2 / (2 * Theta) - mu / 2), 2)) * y
           * sqrt(x2 / 4.0 - y2);
  } else {
    return 0;
  }
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

double HFUtil::Ssf::get() {
  assert(in->getDegeneracy() > 0.0);
  compute(in->getDimension());
  return res;
}

// Compute for 3D systems
void HFUtil::Ssf::compute3D() {
  assert(in->getDegeneracy() > 0.0);
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(yMin, yMax));
  res = 1.0 + itg->getSolution();
}

// Compute for 2D systems
void HFUtil::Ssf::compute2D() {
  const double Theta = in->getDegeneracy();
  assert(in->getDegeneracy() > 0.0);
  auto func1 = [&](const double &y) -> double { return integrand2DOut(y); };
  auto func2 = [&](const double &p) -> double { return integrand2DIn(p); };
  itg2->compute(func1, func2, Itg2DParam(yMin, yMax, 0, M_PI), itgGrid);
  res = itg2->getSolution() + Theta * idr(grid_val, 0);
}

// 3D Integrand
double HFUtil::Ssf::integrand(const double &y) const {
  const double Theta = in->getDegeneracy();
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

// 2D Integrands
double HFUtil::Ssf::integrand2DOut(const double &y) const {
  const double Theta = in->getDegeneracy();
  const double y2 = y * y;
  return 2.0 * y / (exp(y2 / Theta - mu) * M_PI + M_PI);
}

double HFUtil::Ssf::integrand2DIn(const double &p) const {
  const double &Theta = in->getDegeneracy();
  const double y = itg2->getX();
  const double x2 = x * x;
  const double arg = x2 / (2 * Theta) + x * y / Theta * cos(p);
  return SpecialFunctions::coth(arg) - (1.0 / arg);
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