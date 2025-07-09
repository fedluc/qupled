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
using Itg2DParam = Integrator2D::Param;

// // Compute chemical potential
// void HF::computeChemicalPotential2D() {
//   if (in().getDegeneracy() == 0.0) return;
//   ChemicalPotential mu_(in().getDegeneracy());
//   mu_.compute2D();
//   mu = mu_.get();
// }

// Compute ideal density response
void HFUtil::Idr::compute2D() {
  computeIdrFinite2D();
}

void HFUtil::Idr::computeIdrFinite2D() {
  const size_t nx = wvg.size();
  for (size_t i = 0; i < nx; ++i) {
    idr.fill(i, get2D(wvg[i]));
  }
}

// void HF2D::computeIdrGround2D() {
//   const size_t nx = idr.size(0);
//   const size_t nl = idr.size(1);
//   for (size_t i = 0; i < nx; ++i) {
//     for (size_t l = 0; l < nl; ++l) {
//       HFUtil::IdrGround idrTmp(wvg[i], l);
//       idr(i, l) = idrTmp.get();
//     }
//   }
// }

void HF::computeSsf2D() {
  // (in().getDegeneracy() == 0.0) ? computeSsfGround2D() : computeSsfFinite2D();
  computeSsfFinite2D();
}

void HF::computeSsfFinite2D() {
  const bool segregatedItg = in().getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  shared_ptr<Integrator2D> itg2 =
        make_shared<Integrator2D>(in().getIntError());
  for (size_t i = 0; i < wvg.size(); ++i) {
    HFUtil::Ssf2D ssfTmp2D(
        wvg[i], in().getDegeneracy(), mu, wvg.front(), wvg.back(), itgGrid, itg2);
    ssf[i] = ssfTmp2D.get2D() + in().getDegeneracy() * idr(i, 0);
  }
}

// void HF2D::computeSsfGround2D() {
//   for (size_t i = 0; i < wvg.size(); ++i) {
//     HFUtil::SsfGround ssfTmp(wvg[i]);
//     ssf[i] = ssfTmp.get();
//   }
// }

void HF::computeLfc2D() {
  assert(lfc.size() == wvg.size());
  for (auto &s : lfc) {
    s = 1;
  }
}

// Getters
vector<double> HF::getRdf2D(const vector<double> &r) const {
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the radial distribution function");
    return vector<double>();
  }
  return computeRdf2D(r, wvg, ssf);
}

vector<double> HF::getSdr2D() const {
  // need to implement 2D version of Sdr
  const double theta = in().getDegeneracy();
  if (isnan(theta) || theta == 0.0) { return vector<double>(); }
  vector<double> sdr(wvg.size(), -1.5 * theta);
  const double fact = 4 * numUtil::lambda * in().getCoupling() / M_PI;
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x2 = wvg[i] * wvg[i];
    const double phi0 = idr(i, 0);
    sdr[i] *= phi0 / (1.0 + fact / x2 * (1.0 - lfc(i, 0)) * phi0);
  }
  return sdr;
}

double HF::getUInt2D() const {
  if (wvg.size() < 3 || ssf.size() < 3) {
    throwError("No data to compute the internal energy");
    return numUtil::Inf;
  }
  // need to impelement 2D version of computeInternalEnergy
  return computeInternalEnergy(wvg, ssf, in().getCoupling());
}

// -----------------------------------------------------------------
// Idr2D class
// -----------------------------------------------------------------

// Integrand for frequency = l and wave-vector = x
double HFUtil::Idr::integrand2D(const double &x, const double &y, const int &l) const {
  double phi;
  double y2 = y * y;
  double x2 = x * x;
  double x4 = x2 * x2;
  double plT = M_PI * l * Theta;
  double plT2 = plT * plT;
  double exp1 = x4 / 4.0 - x2 * y2 - plT2;
  if (exp1 > 0.0) {
    phi = atan(x2 * plT / exp1)/2.0;
  } else {
    phi = M_PI/2.0 - atan(x2 * plT / exp1)/2.0;
  }
  if (x > 0.0) {
    return y / (exp(y2 / Theta - mu) + 1.0)
           * 2.0 * abs(cos(phi))/ pow((exp1 * exp1 + x4 * plT2), 0.25);
  } else {
    return 0;
  }
}

// Integrand for frequency = 0 and vector = x
double HFUtil::Idr::integrand2D(const double &x, const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  if (x > 0.0) {
    return 1.0 / (Theta * x * pow(cosh(y2 / (2 * Theta) - mu/2), 2)) * y * sqrt(x2 / 4.0 - y2);
  } else {
    return 0; 
  }
}
// Get result of integration
vector<double> HFUtil::Idr::get2D(const double &x) const {
  assert(Theta > 0.0);
  vector<double> res(nl);
  for (int l = 0; l < nl; ++l) {
    auto func = [&](const double &y) -> double {
      return (l == 0) ? integrand2D(x, y) : integrand2D(x, y, l);
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
  return res;
}

// -----------------------------------------------------------------
// IdrGround2D class
// -----------------------------------------------------------------

// Get
// double HFUtil2D::IdrGround2D::get2D() const {
//   const double x2 = x * x;
//   const double Omega2 = Omega * Omega;
//   const double tx = 2.0 * x;
//   const double x2ptx = x2 + tx;
//   const double x2mtx = x2 - tx;
//   const double x2ptx2 = x2ptx * x2ptx;
//   const double x2mtx2 = x2mtx * x2mtx;
//   const double logarg = (x2ptx2 + Omega2) / (x2mtx2 + Omega2);
//   const double part1 = (0.5 - x2 / 8.0 + Omega2 / (8.0 * x2)) * log(logarg);
//   const double part2 =
//       0.5 * Omega * (atan(x2ptx / Omega) - atan(x2mtx / Omega));
//   if (x > 0.0) { return (part1 - part2 + x) / tx; }
//   return 0;
// }

// -----------------------------------------------------------------
// SsfHF2D class
// -----------------------------------------------------------------


inline double coth(double x) {
  return 1.0 / tanh(x);
}

// Integrand
double HFUtil::Ssf2D::integrandOut(const double &y) const {
const double y2 = y * y;
  return 2.0 * y / (exp(y2 / Theta - mu2D) * M_PI + M_PI);
}

double HFUtil::Ssf2D::integrandIn(const double &p) const {
  const double y = itg2->getX();
  const double x2 = x * x;
  const double arg = x2 / (2 * Theta) + x * y / Theta * cos(p);
  return coth(arg) - (1.0 / arg); 
}
// Get result of integration
double HFUtil::Ssf2D::get2D() const {
  assert(Theta > 0.0);
  auto func1 = [&](const double &y) -> double {
    return integrandOut(y);
  };
  auto func2 = [&](const double &p) -> double {
    return integrandIn(p);
  };
  itg2->compute(
      func1,
      func2,
      Itg2DParam(yMin, yMax, 0, M_PI),
      itgGrid);
  return itg2->getSolution();
}

// -----------------------------------------------------------------
// SsfHFGround2D class
// -----------------------------------------------------------------

// Static structure factor at zero temperature
// double HFUtil::SsfGround::get() const {
//   if (x < 2.0) {
//     return (x / 16.0) * (12.0 - x * x);
//   } else {
//     return 1.0;
//   }
// }