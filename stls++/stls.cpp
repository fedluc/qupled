#include <omp.h>
#include "stls.hpp"
#include "chemicalpotential.hpp"

using namespace vecUtil;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

void Stls::compute(){
  buildWvGrid();
  computeChemicalPotential();
  computeIdr();
  computeSsfHF();
  doIterations();
}

// Set up wave-vector grid
void Stls::buildWvGrid(){
  cout << "Assembling wave vector grid: ";
  wvg.push_back(0.0);
  const shared_ptr<StaticInput> &inStat = in.getStaticInput();
  const double dx = inStat->getWaveVectorGridRes();
  const double xmax = inStat->getWaveVectorGridCutoff();
  while(wvg.back() < xmax){
    wvg.push_back(wvg.back() + dx);
  }
  cout << "Done" << endl;
}

// Compute chemical potential
void Stls::computeChemicalPotential(){
  if (in.getDegeneracy() == 0.0) return;
  cout << "Computing chemical potential: "; 
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const vector<double> &guess = statIn->getChemicalPotentialGuess();
  ChemicalPotential mu_(in.getDegeneracy());
  try {
    mu_.compute(guess);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
  mu = mu_.get();
  computedChemicalPotential = true;
  cout << "Done" << endl;
}

// Compute ideal density response
void Stls::computeIdr(){
  if (in.getDegeneracy() == 0.0) return;
  cout << "Computing ideal density response: "; 
  assert(computedChemicalPotential);
  assert(itg != NULL);
  const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>();
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int nx = wvg.size();
  const int nl = statIn->getNMatsubara();
  idr.resize(nx);
  for (int i=0; i<nx; ++i){
    Idr idrTmp(nl, wvg[i], in.getDegeneracy(), mu,
	       wvg.front(), wvg.back(), itg);
    idr[i] = idrTmp.get();
  }
  cout << "Done" << endl;
}

// Compute Hartree-Fock static structure factor
void Stls::computeSsfHF(){
  cout << "Computing HF static structure factor: "; 
  const int nx = wvg.size();
  ssfHF.resize(nx);
  if (in.getDegeneracy() > 0) computeSsfHFFinite();
  else computeSsfHFGround();
  cout << "Done" << endl;
}

void Stls::computeSsfHFFinite(){
  assert(computedChemicalPotential);
  for (int i=0; i<wvg.size(); ++i) {
    SsfHF ssfTmp(wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    ssfHF[i] = ssfTmp.get();
  }
}

void Stls::computeSsfHFGround(){
  for (int i=0; i<wvg.size(); ++i) {
    SsfHF ssfTmp(wvg[i]);
    ssfHF[i] = ssfTmp.get();
  }
}

// Compute static structure factor
void Stls::computeSsf(){
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  if (in.getDegeneracy() > 0) computeSsfFinite();
  else computeSsfGround();
}

// Compute static structure factor at finite temperature
void Stls::computeSsfFinite(){
  assert(computedChemicalPotential);
  assert(slfc.size() > 0);
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], idr[i], slfcOld[i]);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Stls::computeSsfGround(){
  assert(itg != NULL);
  assert(slfc.size() > 0);
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    const double x = wvg[i];
    double yMin = 0.0;
    if (x > 2.0) yMin = x * (x - 2.0);
    const double yMax = x * (x + 2.0);
    Ssf ssfTmp(x, rs, ssfHF[i], slfcOld[i], yMin, yMax, itg);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void Stls::computeSlfc(){
  const int nx = wvg.size();
  assert(ssf.size() == nx);
  assert(slfc.size() == nx);
  const shared_ptr<Interpolator> itp = make_shared<Interpolator>(wvg,ssf);
  for (int i=0; i<nx; ++i){
    Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itg, itp);
    slfc[i] = slfcTmp.get();
  }
}

// stls iterations
void Stls::doIterations() {
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int maxIter = statIn->getNIter();
  const double minErr = statIn->getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  if (verbose) cout << "Structural properties calculation..." << endl;
  while (counter < maxIter && err > minErr ) {
    // Start timing
    double tic = omp_get_wtime();
    // Update SSF
    computeSsf();
    // Update SLFC
    computeSlfc();
    // Update diagnostic
    counter++;
    err = computeError();
    updateSolution();
    // End timing
    double toc = omp_get_wtime();
    // Print diagnostic
    if (verbose) {
       printf("--- iteration %d ---\n", counter);
       printf("Elapsed time: %f seconds\n", toc - tic);
       printf("Residual error: %.5e\n", err);
       fflush(stdout);
     }
  }
  if (verbose) cout << "Done" << endl;
}

// Initial guess for stls iterations
void Stls::initialGuess() {
  const int nx = wvg.size();
  slfcOld.resize(nx);
  slfc.resize(nx);
  fill(slfcOld.begin(), slfcOld.end(), 0.0);
  fill(slfc.begin(), slfc.end(), 1.0);
}

// Compute residual error for the stls iterations
double Stls::computeError(){
  return rms(slfc, slfcOld, false);
}

// Update solution during stls iterations
void Stls::updateSolution(){
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const double aMix = statIn->getMixingParameter();
  slfcOld = sum(mult(slfc, aMix), mult(slfcOld, 1 - aMix));
}

// -----------------------------------------------------------------
// Ideal density response class
// -----------------------------------------------------------------

// Integrand for frequency = l and wave-vector = x
double Idr::integrand(const double y, const int l) {
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
double Idr::integrand(const double y) {
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

// Get at finite temperature
vector<double> Idr::get() {
  assert(Theta > 0.0);
  vector<double> res(nl);
  for (int l=0; l<nl; ++l){
    if (l == 0) {
      auto func = [&](double y)->double{return integrand(y);};
      itg->compute(func, yMin, yMax);
    }
    else {
      auto func = [&](double y)->double{return integrand(y,l);};;
      itg->compute(func, yMin, yMax);
    }
    res[l] = itg->getSolution();
  }
  return res;
}

// Real part at zero temperature
double Idr::re0() const {
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
double Idr::im0() const {
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
double Idr::re0Der() const {
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
// Hartree-Fock static structure factor class
// -----------------------------------------------------------------

// Integrand for finite temperature calculations
double SsfHF::integrand(const double y) {
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

// Get at any temperature
double SsfHF::get() {
  if (Theta == 0.0) return get0();
  auto func = [&](double y)->double{return integrand(y);};
  itg->compute(func, yMin, yMax);
  return 1.0 + itg->getSolution();
}

// Get static structure factor at zero temperature
double SsfHF::get0(){
  assert(Theta == 0.0);
  if (x < 2.0) {
    return (x/16.0)*(12.0 - x*x);
  }
  else {
    return 1.0;
  }
}


// -----------------------------------------------------------------
// Static structure factor class
// -----------------------------------------------------------------

// Get at any temperature
double Ssf::get() {
  if (Theta == 0.0) return get0();
  if (x == 0.0) return 0.0;
  const int nl = idr.size();
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

// Get at zero temperature
double Ssf::get0() {
  if (x == 0.0) return 0.0;
  auto func = [&](double y)->double{return integrand(y);};
  itg->compute(func, yMin, yMax);
  double ssfP;
  try {
    ssfP = plasmon();
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
  return ssfHF + itg->getSolution() + ssfP;
}

// Integrand for zero temperature calculations
double Ssf::integrand(const double Omega) {
  double x2 = x*x;
  double fact = (4.0 * lambda * rs)/(M_PI * x2);
  Idr idrTmp(Omega, x);
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
double Ssf::plasmon() {
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
  RootSolver rsol;
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
double Ssf::drf(const double Omega){
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  const double idrRe = Idr(Omega, x).re0();
  const double wCo = x*x + 2*x;     
  assert(Omega >= wCo);
  return 1.0 + fact * idrRe / (1.0 - fact * slfc * idrRe);
}


// Frequency derivative of the dielectric response function  
double Ssf::drfDer(const double Omega){
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  const Idr idrTmp(Omega, x);
  const double idrRe = idrTmp.re0();
  const double idrReDer = idrTmp.re0Der();
  double denom = (1.0 - fact * slfc * idrRe);
  double w_co = x*x + 2*x;
  assert(Omega >= w_co); 
  return fact * idrReDer / (denom * denom);
}

// -----------------------------------------------------------------
// Static local field correction
// -----------------------------------------------------------------

// Compute static structure factor from interpolator
double Slfc::ssf(const double x_){
  return ssfi->eval(x_);
}

// Get at finite temperature
double Slfc::get() {
  auto func = [&](double y)->double{return integrand(y);};
  itg->compute(func, yMin, yMax);
  return itg->getSolution();
}

// Integrand for finite temperature calculations
double Slfc::integrand(const double y) {
  double y2 = y*y;
  double x2 = x*x;
  if (x > 0.0 && y > 0.0){
    if (x > y){
      return -(3.0/4.0) * y2 * (ssf(y) - 1.0)
	* (1 + (x2 - y2)/(2*x*y)*log((x + y)/(x - y)));
    }
    else if (x < y) {
      return -(3.0/4.0) * y2 * (ssf(y) - 1.0)
	* (1 + (x2 - y2)/(2*x*y)*log((x + y)/(y - x)));
    }
    else {
      return -(3.0/4.0) * y2 * (ssf(y) - 1.0);
    }
  }
  else
    return 0;
}
