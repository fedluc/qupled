#include <omp.h>
#include "stls.hpp"
#include "chemicalpotential.hpp"

using namespace vecUtil;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

void Stls::compute(){
  // Build wave vector grid
  cout << "Assembling wave vector grid: ";
  buildWvGrid();
  cout << "Done." << endl;
  // Chemical potential
  cout << "Computing chemical potential: "; 
  computeChemicalPotential();
  cout << mu << ", Done." << endl;
  // Define integrator object
  itg = make_shared<Integrator1D>();
  // Ideal density response
  cout << "Computing ideal density response: "; 
  computeIdr();
  cout << "Done." << endl;
  // Hartree-Fock static structure factor
  cout << "Computing HF static structure factor: "; 
  computeSsfHF();
  cout << "Done." << endl;
  // Stls iterations
  doIterations();
  // Print results
  for (double el : ssf) cout << el << endl;
}

// Set up wave-vector grid
void Stls::buildWvGrid(){
  wvg.push_back(0.0);
  const shared_ptr<StaticInput> &inStat = in.getStaticInput();
  const double dx = inStat->getWaveVectorGridRes();
  const double xmax = inStat->getWaveVectorGridCutoff();
  while(wvg.back() < xmax){
    wvg.push_back(wvg.back() + dx);
  }
}

// Compute chemical potential
void Stls::computeChemicalPotential(){
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const vector<double> &guess = statIn->getChemicalPotentialGuess();
  ChemicalPotential mu_(in.getDegeneracy());
  mu_.compute(guess);
  mu = mu_.get();
  computedChemicalPotential = true;
}

// Compute ideal density response
void Stls::computeIdr(){
  assert(computedChemicalPotential);
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int nl = statIn->getNMatsubara();
  const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>();
  idr.resize(nl);
  for (int l=0; l<nl; ++l){
    computeIdrSingleFrequency(l);
  }
}

void Stls::computeIdrSingleFrequency(const int l){
  assert(itg != NULL);
  int nx = wvg.size();
  idr[l].resize(nx);
  for (int i=0; i<nx; ++i){
    Idr idrTmp(l, wvg[i], in.getDegeneracy(), mu,
	       wvg.front(), wvg.back(), itg);
    idr[l][i] = idrTmp.get();
  }
}

// Compute Hartree-Fock static structure factor
void Stls::computeSsfHF(){
  assert(computedChemicalPotential);
  const int nx = wvg.size();
  ssfHF.resize(nx);
  for (int i=0; i<nx; ++i){
    SsfHF ssfTmp(wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back(), itg);
    ssfHF[i] = ssfTmp.get();
  }
}

// Compute static structure factor at finite temperature
void Stls::computeSsf(){
  assert(computedChemicalPotential);
  assert(in.getDegeneracy() > 0);
  assert(slfc.size() > 0);
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const double fact1 = 4.0*lambda*rs/M_PI;
  const int nx = wvg.size();
  const int nl = statIn->getNMatsubara();
  // Hartree-Fock (HF) contribution
  ssf = ssfHF;
  // Beyond HF contribution for x = 0
  ssf[0] = 0.0;
  // Beyond HF contribution for x > 0
  for (int i=1; i<nx; ++i){
    const double x2 = wvg[i]*wvg[i];
    double fact2 = 0.0;
    for (int l=0; l<nl; ++l) {
      const double fact3 = 1.0 + fact1/x2*(1- slfc[i])*idr[l][i];
      double fact4 = idr[l][i]*idr[l][i]/fact3;
      if (l>0) fact4 *= 2;
      fact2 += fact4;
    }
    ssf[i] += -1.5*fact1/x2*Theta*(1 - slfc[i])*fact2;
  }
}

// Compute static local field correction
void Stls::computeSlfc(){
  const int nx = wvg.size();
  assert(computedChemicalPotential);
  assert(ssf.size() == nx);
  const shared_ptr<Interpolator> itp = make_shared<Interpolator>(wvg,ssf);
  if (slfc.size() == 0) slfc.resize(nx);
  for (int i=0; i<nx; ++i){
    Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itg, itp);
    slfc[i] = slfcTmp.get();
  }
}

// stls iterations
void Stls::doIterations() {
  const bool verbose = true; // FIX THIS!
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

// Integrand for ideal density response
double Idr::integrand(const double y){
  if (l == 0)
    return x0(y);
  else
    return xl(y);
}


// Integrand for frequency = l and wave-vector = x
double Idr::xl(const double y) {
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
double Idr::x0(const double y) {
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
double Idr::get() {
  auto func = [&](double y)->double{return integrand(y);};
  itg->compute(func, yMin, yMax);
  return itg->getSolution();
}

// Real part at zero temperature
double Idr::re0() {
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
double Idr::im0() {
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
double Idr::re0Der() {
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

// Get at finite temperature
double SsfHF::get() {
  auto func = [&](double y)->double{return integrand(y);};
  itg->compute(func, yMin, yMax);
  return 1.0 + itg->getSolution();
}

// Get static structure factor at zero temperature
double SsfHF::get0(){
  if (x < 2.0) {
    return (x/16.0)*(12.0 - x*x);
  }
  else {
    return 1.0;
  }
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
