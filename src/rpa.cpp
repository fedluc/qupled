#include <omp.h>
#include "util.hpp"
#include "rpa.hpp"
#include "chemicalpotential.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace thermoUtil;
using namespace binUtil;

// Constructor
Rpa::Rpa(const RpaInput &in_,
	 const bool verbose_,
	 const bool doFullCompute) : in(in_),
				     verbose(verbose_),
				     itg(in_.getIntError()) {
  if (verbose) cout << "Assembling wave vector grid: ";
  buildWvGrid();
  if (verbose) cout << "Done" << endl;
  if (verbose) cout << "Computing chemical potential: "; 
  computeChemicalPotential();
  if (verbose) cout << "Done" << endl;
  if (verbose) cout << "Computing ideal density response: "; 
  computeIdr();
  if (verbose) cout << "Done" << endl;
  if (verbose) cout << "Computing HF static structure factor: "; 
  computeSsfHF();
  if (verbose) cout << "Done" << endl;
  if (doFullCompute) {
    if (verbose) cout << "Computing static local field correction: "; 
    computeSlfc();
    if (verbose) cout << "Done" << endl;
    if (verbose) cout << "Computing static structure factor: "; 
    computeSsf();
    if (verbose) cout << "Done" << endl;
  }
}

// Set up wave-vector grid
void Rpa::buildWvGrid(){
  wvg.push_back(0.0);
  const double dx = in.getWaveVectorGridRes();
  const double xmax = in.getWaveVectorGridCutoff();
  while(wvg.back() < xmax){
    wvg.push_back(wvg.back() + dx);
  }
}

// Compute chemical potential
void Rpa::computeChemicalPotential(){
  if (in.getDegeneracy() == 0.0) return;
  const vector<double> &guess = in.getChemicalPotentialGuess();
  ChemicalPotential mu_(in.getDegeneracy());
  try {
    mu_.compute(guess);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
  mu = mu_.get();
}

// Compute ideal density response
void Rpa::computeIdr(){
  if (in.getDegeneracy() == 0.0) return;
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  idr.resize(nx, nl);
  for (int i=0; i<nx; ++i){
    Idr idrTmp(nl, wvg[i], in.getDegeneracy(), mu,
	       wvg.front(), wvg.back(), itg);
    idr.fill(i, idrTmp.get());
  }
}

// Compute Hartree-Fock static structure factor
void Rpa::computeSsfHF(){
  const int nx = wvg.size();
  ssfHF.resize(nx);
  if (in.getDegeneracy() == 0.0) {
    computeSsfHFGround();
    return;
  }
  computeSsfHFFinite();
}

void Rpa::computeSsfHFFinite(){
  for (size_t i=0; i<wvg.size(); ++i) {
    SsfHF ssfTmp(wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    ssfHF[i] = ssfTmp.get();
  }
}

void Rpa::computeSsfHFGround(){
  for (size_t i=0; i<wvg.size(); ++i) {
    SsfHFGround ssfTmp(wvg[i]);
    ssfHF[i] = ssfTmp.get();
  }
}

// Compute static structure factor
void Rpa::computeSsf(){
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  if (in.getDegeneracy() == 0.0) {
    computeSsfGround();
    return;
  }
  computeSsfFinite();
}

// Compute static structure factor at finite temperature
void Rpa::computeSsfFinite(){
  assert(slfc.size() > 0);
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  const int nl = idr.size(1);
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], slfc[i], nl, &idr(i));
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Rpa::computeSsfGround(){
  assert(slfc.size() > 0);
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    const double x = wvg[i];
    double yMin = 0.0;
    if (x > 2.0) yMin = x * (x - 2.0);
    const double yMax = x * (x + 2.0);
    SsfGround ssfTmp(x, rs, ssfHF[i], slfc[i], yMin, yMax, itg);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void Rpa::computeSlfc(){
  slfc.resize(wvg.size());
  for (auto& s : slfc) { s = 0; }
}

// Getters
vector<double> Rpa::getRdf(const vector<double> &r) const {
  return computeRdf(r, wvg, ssf);
}

vector<double> Rpa::getSdr() const {
  if (in.getDegeneracy() == 0.0) {
    std::cout << "The static density response cannot be computed in the ground state." << std::endl;
    return vector<double>();
  }
  vector<double> sdr(wvg.size(), -1.5 * in.getDegeneracy());
  const double fact = 4 *lambda * in.getCoupling() / M_PI;
  for (size_t i=0; i<wvg.size(); ++i){
    const double x2 = wvg[i] * wvg[i];
    const double phi0 = idr(i,0);
    sdr[i] *= phi0/ (1.0 + fact/x2 * (1.0 - slfc[i]) * phi0);
  }
  return sdr;
}

double Rpa::getUInt() const {
  return computeInternalEnergy(wvg, ssf, in.getCoupling());
};  

// -----------------------------------------------------------------
// Idr class
// -----------------------------------------------------------------

// Integrand for frequency = l and wave-vector = x
double Idr::integrand(const double y, const int l) const {
  double y2 = y*y;
  double x2 = x*x;
  double txy = 2*x*y; 
  double tplT = 2*M_PI*l*Theta;
  double tplT2 = tplT*tplT;
  if (x > 0.0) {
    return 1.0/(2*x)*y/(exp(y2/Theta - mu) + 1.0)
      *log(((x2+txy)*(x2+txy) + tplT2)/((x2-txy)*(x2-txy) + tplT2));
  }
  else {
    return 0;
  }
}

// Integrand for frequency = 0 and vector = x
double Idr::integrand(const double y) const {
  double y2 = y*y;
  double x2 = x*x;
  double xy = x*y;
  if (x > 0.0){
    if (x < 2*y){
      return 1.0/(Theta*x)*((y2 - x2/4.0)*log((2*y + x)/(2*y - x)) + xy)
        *y/(exp(y2/Theta - mu) + exp(-y2/Theta + mu) + 2.0);
    }
    else if (x > 2*y){
      return 1.0/(Theta*x)*((y2 - x2/4.0)*log((2*y + x)/(x  - 2*y)) + xy)
        *y/(exp(y2/Theta - mu) + exp(-y2/Theta + mu) + 2.0);
    }
    else {
      return 1.0/(Theta)*y2/(exp(y2/Theta - mu) + exp(-y2/Theta + mu) + 2.0);;
    }
  }
  else{
    return (2.0/Theta)*y2/(exp(y2/Theta - mu) + exp(-y2/Theta + mu) + 2.0);
  }
}

// Get result of integration
vector<double> Idr::get() const {
  assert(Theta > 0.0);
  vector<double> res(nl);
  for (int l=0; l<nl; ++l){
    if (l == 0) {
      auto func = [&](double y)->double{return integrand(y);};
      itg.compute(func, yMin, yMax);
    }
    else {
      auto func = [&](double y)->double{return integrand(y,l);};;
      itg.compute(func, yMin, yMax);
    }
    res[l] = itg.getSolution();
  }
  return res;
}

// -----------------------------------------------------------------
// IdrGround class
// -----------------------------------------------------------------

// Real part at zero temperature
double IdrGround::re0() const {
  double adder1 = 0.0;
  double adder2 = 0.0;
  double preFactor = 0.0;
  if (x > 0.0) {
    double x_2 = x/2.0;
    double Omega_2x = Omega/(2.0*x);
    double sumFactor = x_2 + Omega_2x;
    double diffFactor = x_2 - Omega_2x;
    double sumFactor2 = sumFactor*sumFactor;
    double diffFactor2 = diffFactor*diffFactor;
    preFactor = 0.5;
    if (sumFactor != 1.0) {
      double log_sum_arg = (sumFactor + 1.0)/(sumFactor - 1.0);
      if (log_sum_arg < 0.0) log_sum_arg = -log_sum_arg;
      adder1 = 1.0/(4.0*x)*(1.0 - sumFactor2)*log(log_sum_arg);
    }
    if (diffFactor != 1.0 && diffFactor != -1.0) {
      double log_diff_arg = (diffFactor + 1.0)/(diffFactor - 1.0);
      if (log_diff_arg < 0.0) log_diff_arg = -log_diff_arg;
      adder2 = 1.0/(4.0*x)*(1.0 - diffFactor2)*log(log_diff_arg);
    }
  }
  return preFactor + adder1 + adder2;
}

// Imaginary part at zero temperature
double IdrGround::im0() const {
  double preFactor = 0.0;
  double adder1 = 0.0;
  double adder2 = 0.0;
  if (x > 0.0) {
    double x_2 = x/2.0;
    double Omega_2x = Omega/(2.0*x);
    double sumFactor = x_2 + Omega_2x;
    double diffFactor = x_2 - Omega_2x;
    double sumFactor2 = sumFactor*sumFactor;
    double diffFactor2 = diffFactor*diffFactor;
    preFactor = -M_PI/(4.0*x);
    if (sumFactor2 < 1.0) {
      adder1 = 1 - sumFactor2;
    }
    if (diffFactor2 < 1.0) {
      adder2 = 1 - diffFactor2;
    }
  }
  return preFactor * (adder1 - adder2);
}

// Frequency derivative of the real part at zero temperature
double IdrGround::re0Der() const {
  double adder1 = 0.0;
  double adder2 = 0.0;
  double x_2 = x/2.0;
  double Omega_2x  = Omega/(2.0*x);
  double sumFactor = x_2 + Omega_2x;
  double diffFactor = x_2 - Omega_2x;
  if (sumFactor != 1.0) {
    double log_sum_arg = (sumFactor + 1.0)/(sumFactor - 1.0);
    if (log_sum_arg < 0.0) log_sum_arg = -log_sum_arg;
    adder1 = 1.0/(4.0*x*x)*(1.0 - sumFactor*log(log_sum_arg));
  }
  if (diffFactor != 1.0 && diffFactor != -1.0) {
    double log_diff_arg = (diffFactor + 1.0)/(diffFactor - 1.0);
    if (log_diff_arg < 0.0) log_diff_arg = -log_diff_arg;
    adder2 = -1.0/(4.0*x*x)*(1.0 - diffFactor*log(log_diff_arg));
  }
  return adder1 + adder2;
}


// -----------------------------------------------------------------
// SsfHF class
// -----------------------------------------------------------------

// Integrand
double SsfHF::integrand(const double y) const {
  double y2 = y*y;
  double ypx = y + x;
  double ymx = y - x;
  if (x > 0.0){
    return -3.0*Theta/(4.0*x)*y/(exp(y2/Theta - mu) + 1.0)
      *log((1 + exp(mu - ymx*ymx/Theta))/(1 + exp(mu - ypx*ypx/Theta)));
  }
  else {
    return -3.0*y2/((1.0 + exp(y2/Theta - mu))*(1.0 + exp(y2/Theta - mu)));
  }
}

// Get result of integration
double SsfHF::get() const {
  assert(Theta > 0.0);
  auto func = [&](double y)->double{return integrand(y);};
  itg.compute(func, yMin, yMax);
  return 1.0 + itg.getSolution();
}

// -----------------------------------------------------------------
// SsfHFGround class
// -----------------------------------------------------------------

// Static structure factor at zero temperature
double SsfHFGround::get() const {
  if (x < 2.0) {
    return (x/16.0)*(12.0 - x*x);
  }
  else {
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
  const double fact1 = 4.0*lambda*rs/M_PI;
  const double x2 = x*x;
  double fact2 = 0.0;
  for (int l=0; l<nl; ++l) {
    const double fact3 = 1.0 + fact1/x2*(1- slfc)*idr[l];
    double fact4 = idr[l]*idr[l]/fact3;
    if (l>0) fact4 *= 2;
    fact2 += fact4;
  }
  return ssfHF - 1.5 * fact1/x2 * Theta * (1 - slfc) * fact2;
}

// -----------------------------------------------------------------
// SsfGround class
// -----------------------------------------------------------------

// Get result of integration
double SsfGround::get() const {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](double y)->double{return integrand(y);};
  itg.compute(func, yMin, yMax);
  double ssfP;
  try {
    ssfP = plasmon();
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
  return ssfHF + itg.getSolution() + ssfP;
}

// Integrand for zero temperature calculations
double SsfGround::integrand(const double Omega) const {
  double x2 = x*x;
  double fact = (4.0 * lambda * rs)/(M_PI * x2);
  IdrGround idrTmp(Omega, x);
  const double idrRe = idrTmp.re0();
  const double idrIm = idrTmp.im0();
  const double factRe = 1 + fact * (1 - slfc) * idrRe;
  const double factIm = fact * (1 - slfc) * idrIm;
  const double factRe2 = factRe * factRe;
  const double factIm2 = factIm * factIm;
  return 1.5/(M_PI)* idrIm * (1.0/(factRe2 + factIm2) - 1.0);
}

// NOTE: At the plasmon frequency, the imaginary part of the ideal
// density response is zero. Hence, all the following definitions
// for the dielectric function are constructed with the assumption
// that the imaginary part of the ideal density response is zero
// and should not be used for frequencies omega < 2*xx + xx^2 (with
// xx a normalized wave-vector) where such approximation is not
// valid

// Plasmon contribution to the static structure factor
double SsfGround::plasmon() const {
  // Look for a region where the dielectric function changes sign
  bool search_root = false;
  const double wCo = x*x + 2*x;
  const double dw = wCo;
  const double wLo = wCo;
  double wHi;
  const int signLo = (drf(wLo) >= 0) ? 1 : -1;
  for (int i=1; i<1000; i++) {
    wHi = wLo + dw * i;
    const double signHi = (drf(wHi) >= 0) ? 1 : -1;;
    if (signHi != signLo) {
      search_root = true;
      break;
    }
  }
  // Return if no root can be found
  if (!search_root) return 0;
  // Compute plasmon frequency
  auto func = [this](double Omega)->double{return drf(Omega);};
  const double guess[] = {wLo, wHi};
  BrentRootSolver rsol;
  rsol.solve(func, vector<double>(begin(guess),end(guess)));
  if (!rsol.success()) {
    throw runtime_error("Plasmon solver: the root solver "
			"did not converge to the desired accuracy.");
  }
  // Output
  const double fact = (4.0 *lambda *rs)/(M_PI * x * x);
  return 1.5 / (fact * abs(drfDer(rsol.getSolution())));
}

// Dielectric response function
double SsfGround::drf(const double Omega) const {
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  const double idrRe = IdrGround(Omega, x).re0();     
  assert(Omega >= x*x + 2*x);
  return 1.0 + fact * idrRe / (1.0 - fact * slfc * idrRe);
}


// Frequency derivative of the dielectric response function  
double SsfGround::drfDer(const double Omega) const {
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  Integrator1D itgTmp = itg;
  const IdrGround idrTmp(Omega, x);
  const double idrRe = idrTmp.re0();
  const double idrReDer = idrTmp.re0Der();
  double denom = (1.0 - fact * slfc * idrRe);
  assert(Omega >= x*x + 2*x); 
  return fact * idrReDer / (denom * denom);
}

