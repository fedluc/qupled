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
Rpa::Rpa(const std::shared_ptr<const Input> &in_, const bool verbose_)
    : HF(in_, verbose_) {
  // Allocate arrays to the correct size
  const size_t nx = wvg.size();
  const size_t nl = in().getNMatsubara();
  idr.resize(nx, nl);
  ssfHF.resize(nx);
}

// Initialize basic properties
void Rpa::init() {
  HF::init();
  print("Computing Hartree-Fock static structure factor: ");
  computeSsfHF();
  println("Done");
}

// Compute Hartree-Fock static structure factor
void Rpa::computeSsfHF() {
  HF hf(inPtr, false);
  hf.compute();
  ssfHF = hf.getSsf();
}

// Compute static structure factor at finite temperature
void Rpa::computeSsfFinite() {
  const size_t nx = wvg.size();
  for (size_t i = 0; i < nx; ++i) {
    RpaUtil::Ssf ssfTmp(wvg[i], ssfHF[i], lfc[i], inPtr, idr[i]);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Rpa::computeSsfGround() {
  const size_t nx = wvg.size();
  for (size_t i = 0; i < nx; ++i) {
    const double x = wvg[i];
    RpaUtil::SsfGround ssfTmp(x, ssfHF[i], lfc[i], itg, inPtr);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void Rpa::computeLfc() {
  assert(lfc.size() == wvg.size());
  for (auto &s : lfc) {
    s = 0;
  }
}

// -----------------------------------------------------------------
// SsfBase class
// -----------------------------------------------------------------

// Normalized interaction potential
double RpaUtil::SsfBase::ip() const {
  const double rs = in->getCoupling();
  if (in->getDimension() == dimensionsUtil::Dimension::D2) {
    return sqrt(2.0) * rs / x;
  } else {
    return 4.0 * numUtil::lambda * rs / (M_PI * x * x);
  }
}

// -----------------------------------------------------------------
// Ssf class
// -----------------------------------------------------------------

double RpaUtil::Ssf::get() {
  assert(in->getDegeneracy() > 0.0);
  if (x == 0.0) return 0.0;
  if (in->getCoupling() == 0.0) return ssfHF;
  compute(in->getDimension());
  return res;
}

void RpaUtil::Ssf::compute3D() {
  const double Theta = in->getDegeneracy();
  const double isStatic = lfc.size() == 1;
  double suml = 0.0;
  for (size_t l = 0; l < idr.size(); ++l) {
    const double &idrl = idr[l];
    const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
    const double denom = 1.0 + ip() * idrl * (1 - lfcl);
    const double f = idrl * idrl * (1 - lfcl) / denom;
    suml += (l == 0) ? f : 2 * f;
  }
  res = ssfHF - 1.5 * ip() * Theta * suml;
}

void RpaUtil::Ssf::compute2D() {
  const double Theta = in->getDegeneracy();
  const double isStatic = lfc.size() == 1;
  double suml = 0.0;
  for (size_t l = 0; l < idr.size(); ++l) {
    const double &idrl = idr[l];
    const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
    const double denom = 1.0 + ip() * idrl * (1 - lfcl);
    const double f = idrl * idrl * (1 - lfcl) / denom;
    suml += (l == 0) ? f : 2 * f;
  }
  res = ssfHF - ip() * Theta * suml;
}

// -----------------------------------------------------------------
// SsfGround class
// -----------------------------------------------------------------

double RpaUtil::SsfGround::get() {
  const double OmegaMax = in->getFrequencyCutoff();
  const double rs = in->getCoupling();
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(0, OmegaMax));
  return 1.5 / (M_PI)*itg->getSolution() + ssfHF;
}

double RpaUtil::SsfGround::integrand(const double &Omega) const {
  const double idr = HFUtil::IdrGround(x, Omega).get();
  return idr / (1.0 + ip() * idr * (1.0 - lfc[0])) - idr;
}
