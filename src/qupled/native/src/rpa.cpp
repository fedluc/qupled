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
Rpa::Rpa(const Input &in_, const bool verbose_)
    : HF(in_, verbose_) {
  // Allocate arrays to the correct size
  const size_t nx = wvg.size();
  const size_t nl = in.getNMatsubara();
  idr.resize(nx, nl);
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
  HF::init();
  print("Computing Hartree-Fock static structure factor: ");
  computeSsfHF();
  println("Done");
}

// Compute Hartree-Fock static structure factor
void Rpa::computeSsfHF() {
  HF hf(in, false);
  hf.compute();
  ssfHF = hf.getSsf();
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
  assert(slfc.size() == nx);
  assert(ssf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    RpaUtil::Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], slfc[i], idr[i]);
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
    RpaUtil::SsfGround ssfTmp(x, rs, ssfHF[i], slfc[i], OmegaMax, itg);
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

// -----------------------------------------------------------------
// Ssf class
// -----------------------------------------------------------------

// Get at finite temperature for any scheme
double RpaUtil::Ssf::get() const {
  assert(Theta > 0.0);
  if (rs == 0.0) return ssfHF;
  if (x == 0.0) return 0.0;
  double fact2 = 0.0;
  for (size_t l = 0; l < idr.size(); ++l) {
    const double &idrl = idr[l];
    const double fact3 = 1.0 + ip * (1 - slfc) * idrl;
    double fact4 = idrl * idrl / fact3;
    if (l > 0) fact4 *= 2;
    fact2 += fact4;
  }
  return ssfHF - 1.5 * ip * Theta * (1 - slfc) * fact2;
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
  return idr / (1.0 + ip * idr * (1.0 - slfc)) - idr;
}
