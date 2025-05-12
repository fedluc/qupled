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
  const double Theta = in().getDegeneracy();
  const double rs = in().getCoupling();
  const size_t nx = wvg.size();
  for (size_t i = 0; i < nx; ++i) {
    RpaUtil::Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], lfc[i], idr[i]);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Rpa::computeSsfGround() {
  const double rs = in().getCoupling();
  const double OmegaMax = in().getFrequencyCutoff();
  const size_t nx = wvg.size();
  for (size_t i = 0; i < nx; ++i) {
    const double x = wvg[i];
    RpaUtil::SsfGround ssfTmp(x, rs, ssfHF[i], lfc[i], OmegaMax, itg);
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
// Ssf class
// -----------------------------------------------------------------

// Get at finite temperature for any scheme
double RpaUtil::Ssf::get() const {
  if (rs == 0.0) return ssfHF;
  if (x == 0.0) return 0.0;
  const double isStatic = lfc.size() == 1;
  double suml = 0.0;
  for (size_t l = 0; l < idr.size(); ++l) {
    const double &idrl = idr[l];
    const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
    const double denom = 1.0 + ip * idrl * (1 - lfcl);
    const double f = idrl * idrl * (1 - lfcl) / denom;
    suml += (l == 0) ? f : 2 * f;
  }
  return ssfHF - 1.5 * ip * Theta * suml;
}

// -----------------------------------------------------------------
// SsfGround class
// -----------------------------------------------------------------

double RpaUtil::SsfGround::get() {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(0, OmegaMax));
  return 1.5 / (M_PI)*itg->getSolution() + ssfHF;
}

double RpaUtil::SsfGround::integrand(const double &Omega) const {
  const double idr = HFUtil::IdrGround(x, Omega).get();
  return idr / (1.0 + ip * idr * (1.0 - lfc[0])) - idr;
}
