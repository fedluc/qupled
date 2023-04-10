#include <omp.h>
#include "util.hpp"
#include "stls.hpp"
#include "chemicalpotential.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace thermoUtil;
using namespace binUtil;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

void Stls::compute(){
  init();  
  if (verbose) cout << "Structural properties calculation ..." << endl;
  doIterations();
  if (verbose) cout << "Done" << endl;
}

// Initialization
void Stls::init(){
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
}

// Set up wave-vector grid
void Stls::buildWvGrid(){
  wvg.push_back(0.0);
  const double dx = in.getWaveVectorGridRes();
  const double xmax = in.getWaveVectorGridCutoff();
  while(wvg.back() < xmax){
    wvg.push_back(wvg.back() + dx);
  }
}

// Compute chemical potential
void Stls::computeChemicalPotential(){
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
  computedChemicalPotential = true;
}

// Compute ideal density response
void Stls::computeIdr(){
  if (in.getDegeneracy() == 0.0) return;
  assert(computedChemicalPotential);
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
void Stls::computeSsfHF(){
  const int nx = wvg.size();
  ssfHF.resize(nx);
  if (in.getDegeneracy() > 0) computeSsfHFFinite();
  else computeSsfHFGround();
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
    SsfHFGround ssfTmp(wvg[i]);
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
  const int nl = idr.size(1);
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], slfcOld[i], nl, &idr(i));
    ssf[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Stls::computeSsfGround(){
  assert(slfc.size() > 0);
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    const double x = wvg[i];
    double yMin = 0.0;
    if (x > 2.0) yMin = x * (x - 2.0);
    const double yMax = x * (x + 2.0);
    SsfGround ssfTmp(x, rs, ssfHF[i], slfcOld[i], yMin, yMax, itg);
    ssf[i] = ssfTmp.get();
  }
}

// Compute static local field correction
void Stls::computeSlfc(){
  const int nx = wvg.size();
  assert(ssf.size() == nx);
  assert(slfc.size() == nx);
  computeSlfcStls();
  if (useIet) computeSlfcIet();
}

void Stls::computeSlfcStls() {
  const int nx = wvg.size();
  const Interpolator1D itp(wvg,ssf);
  for (int i=0; i<nx; ++i) {
    Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itp, itg);
    slfc[i] = slfcTmp.get();
  }
}

void Stls::computeSlfcIet() {
   Integrator2D itg2;
   const bool segregatedItg = in.getInt2DScheme() == "segregated";
   const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
   const Interpolator1D ssfItp(wvg, ssf);
   const Interpolator1D slfcItp(wvg, slfcOld);
   if (bf.size() == 0) computeBf();
   const Interpolator1D bfItp(wvg, bf);
   for (int i=0; i<wvg.size(); ++i){
     SlfcIet slfcTmp(wvg[i], wvg.front(), wvg.back(),
		     ssfItp, slfcItp, bfItp, itgGrid, itg2);
     slfc[i] += slfcTmp.get();
   }
}

// Compute bridge function
void Stls::computeBf() {
  const int nx = wvg.size();
  Integrator1DFourier itgF(0, 1e-10);
  bf.resize(nx);
  for (int i=0; i<nx; ++i){ 
    BridgeFunction bfTmp(in.getTheory(), in.getIETMapping(),
			 in.getCoupling(), in.getDegeneracy(),
			 wvg[i], itgF);
    bf[i] = bfTmp.get();
  }
}

// stls iterations
void Stls::doIterations() {
  const int maxIter = in.getNIter();
  const int outIter = in.getOutIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  while (counter < maxIter+1 && err > minErr ) {
    // Start timing
    double tic = omp_get_wtime();
    // Update static structure factor
    computeSsf();
    // Update static local field correction
    computeSlfc();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // Write output
    if (counter % outIter == 0) { writeOutput();};
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
}

// Initial guess for stls iterations
void Stls::initialGuess() {
  const int nx = wvg.size();
  slfcOld.resize(nx);
  slfc.resize(nx);
  // From recovery file
  if (in.getRestartFileName() != EMPTY_STRING) {
    vector<double> wvgFile;
    vector<double> slfcFile;
    readRestart(wvgFile, slfcFile);
    const Interpolator1D slfci(wvgFile, slfcFile);
    const double xmaxi = wvgFile.back();
    for (int i=0; i<wvg.size(); ++i) {
      const double x = wvg[i];
      if (x <= xmaxi) { slfcOld[i] = slfci.eval(x);}
      else { slfcOld[i] = 1.0; }
    }
    return;
  }
  // From guess in input
  if (in.getGuess().wvg.size() > 0) {
    const Interpolator1D slfci(in.getGuess().wvg, in.getGuess().property);
    const double xmaxi = in.getGuess().wvg.back();
    for (int i=0; i<wvg.size(); ++i) {
      const double x = wvg[i];
      if (x <= xmaxi) { slfcOld[i] = slfci.eval(x);}
      else { slfcOld[i] = 1.0; }
    }
    return;
  }  
  // Default
  fill(slfcOld.begin(), slfcOld.end(), 0.0);
}

// Compute residual error for the stls iterations
double Stls::computeError(){
  return rms(slfc, slfcOld, false);
}

// Update solution during stls iterations
void Stls::updateSolution(){
  const double aMix = in.getMixingParameter();
  slfcOld = sum(mult(slfc, aMix), mult(slfcOld, 1 - aMix));
}

// Getters
vector<double> Stls::getRdf(const vector<double> &r) const {
  assert(ssf.size() > 0 && wvg.size() > 0);
  const Interpolator1D itp(wvg, ssf);
  const int nr = r.size();
  vector<double> rdf(nr);
  Integrator1DFourier itgF(0.0);
  for (int i=0; i<nr; ++i){
    const Rdf rdfTmp(r[i], wvg.back(), itp, itgF);
    rdf[i] = rdfTmp.get();
  }
  return rdf;
}

vector<double> Stls::getSdr() const {
  if (in.getDegeneracy() == 0.0) {
    throw runtime_error("The static density response cannot be computed in the ground state.");
  };
  vector<double> sdr(wvg.size());
  const double fact = 4 *lambda * in.getCoupling() / M_PI;
  const double Theta = -1.5 * in.getDegeneracy();
  for (int i=0; i<wvg.size(); ++i){
    sdr[i] = idr(i,0)/ (1.0 + fact/(wvg[i] * wvg[i]) * (1.0 - slfc[i]) * idr(i,0));
  }
  transform(sdr.begin(), sdr.end(), sdr.begin(), [&Theta](double el){return Theta*el;});
  return sdr;
}

double Stls::getUInt() const {
  const Interpolator1D itp(wvg, ssf);
  Integrator1D itgTmp;
  const InternalEnergy uInt(in.getCoupling(), wvg.front(), wvg.back(), itp, itgTmp);
  return uInt.get();
}

// Write output files
void Stls::writeOutput() const{
  if (!writeFiles) return; 
  writeRestart();
}

// Restart files
void Stls::writeRestart() const {
  const string fileName = format<double,double>("restart_rs%.3f_theta%.3f_"
						+ in.getTheory() + ".bin",
						in.getCoupling(),
						in.getDegeneracy());
  ofstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Output file " + fileName + " could not be created.");
  }
  int nx = wvg.size();
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<decltype(wvg)>(file, wvg);
  writeDataToBinary<decltype(slfc)>(file, slfc);
  file.close();
  if (!file) {
    throw runtime_error("Error in writing to file " + fileName);
  }
}

void Stls::readRestart(vector<double> &wvgFile,
		       vector<double> &slfcFile) const {
  const string fileName = in.getRestartFileName();
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Output file " + fileName + " could not be opened.");
  }
  int nx;
  readDataFromBinary<int>(file, nx);
  wvgFile.resize(nx);
  slfcFile.resize(nx);
  readDataFromBinary<decltype(wvgFile)>(file, wvgFile);
  readDataFromBinary<decltype(slfcFile)>(file, slfcFile);
  file.close();
  if (!file) {
    throw runtime_error("Error in reading from file " + fileName);
  }
}

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
double SsfGround::drf(const double Omega) const {
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  const double idrRe = IdrGround(Omega, x).re0();
  const double wCo = x*x + 2*x;     
  assert(Omega >= wCo);
  return 1.0 + fact * idrRe / (1.0 - fact * slfc * idrRe);
}


// Frequency derivative of the dielectric response function  
double SsfGround::drfDer(const double Omega) const {
  const double fact = (4.0 * lambda * rs)/(M_PI * x * x);
  Integrator1D itgTmp;
  const IdrGround idrTmp(Omega, x);
  const double idrRe = idrTmp.re0();
  const double idrReDer = idrTmp.re0Der();
  double denom = (1.0 - fact * slfc * idrRe);
  double w_co = x*x + 2*x;
  assert(Omega >= w_co); 
  return fact * idrReDer / (denom * denom);
}

// -----------------------------------------------------------------
// SlfcBase class
// -----------------------------------------------------------------

// Compute static structure factor from interpolator
double SlfcBase::ssf(const double x_) const {
  return ssfi.eval(x_);
}

// -----------------------------------------------------------------
// Slfc class
// -----------------------------------------------------------------

// Get result of integration
double Slfc::get() const {
  auto func = [&](double y)->double{return integrand(y);};
  itg.compute(func, yMin, yMax);
  return itg.getSolution();
}

// Integrand
double Slfc::integrand(const double y) const {
  double y2 = y*y;
  double x2 = x*x;
  if (x == 0.0 || y == 0.0) { return 0.0; }
  if (x == y) { return -(3.0/4.0) * y2 * (ssf(y) - 1.0); };
  if (x > y){
    return -(3.0/4.0) * y2 * (ssf(y) - 1.0)
      * (1 + (x2 - y2)/(2*x*y)*log((x + y)/(x - y)));
  } 
  return -(3.0/4.0) * y2 * (ssf(y) - 1.0)
    * (1 + (x2 - y2)/(2*x*y)*log((x + y)/(y - x)));
}

// -----------------------------------------------------------------
// SlfcIet class
// -----------------------------------------------------------------

// Compute static local field correction from interpolator
double SlfcIet::slfc(const double x_) const{
  return slfci.eval(x_);
}

// Compute bridge function from interpolator
double SlfcIet::bf(const double x_) const {
  return bfi.eval(x_);
}

// Get at finite temperature
double SlfcIet::get() const {
  if (x == 0.0) return 0.0;
  auto wMin = [&](double y)->double{return (y > x) ? y - x : x - y;};
  auto wMax = [&](double y)->double{return min(yMax, x + y);};
  auto func1 = [&](double y)->double{return integrand1(y);};
  auto func2 = [&](double w)->double{return integrand2(w);};
  itg.compute(func1, func2, yMin, yMax, wMin, wMax, itgGrid);
  return 3.0/(8.0*x) * itg.getSolution() + bf(x);
}

// Level 1 integrand
double SlfcIet::integrand1(const double y) const {
  if (y == 0.0) return 0.0;
  return (-bf(y) - (ssf(y) - 1.0)*(slfc(y) - 1.0)) / y;
}

// Level 2 integrand
double SlfcIet::integrand2(const double w) const {
  const double y = itg.getX();
  const double y2 = y*y;
  const double w2 = w*w;
  const double x2 = x*x;
  return (w2 - y2 - x2) * w * (ssf(w) - 1.0);
}


// -----------------------------------------------------------------
// BridgeFunction class
// -----------------------------------------------------------------

double BridgeFunction::get() const {
  if (theory == "STLS-HNC" || theory == "QSTLS-HNC") { return hnc(); }
  if (theory == "STLS-IOI" || theory == "QSTLS-IOI") { return ioi(); }
  if (theory == "STLS-LCT" || theory == "QSTLS-LCT") { return lct(); }
  throw runtime_error("Unknown theory to compute the bridge function term");
}

double BridgeFunction::couplingParameter() const {
  const double fact = 2 * lambda * lambda * rs;
  if (mapping == "sqrt") { return fact/sqrt(1 + Theta * Theta); }
  if (mapping == "linear") { return fact/(1 + Theta); }
  if (Theta != 0.0) { return fact/Theta; }
  throw runtime_error("The standard iet mapping cannot be used in the "
		      "ground state");
}

double BridgeFunction::hnc() const {
  return 0.0;
}

double BridgeFunction::ioi() const {
  const double l2 = lambda*lambda;
  const double l3 = l2*lambda;
  const double l4 = l3*lambda;
  const double l5 = l4*lambda;
  const double l6 = l5*lambda;
  const double l7 = l6*lambda;
  const double l8 = l7*lambda;
  const double Gamma = couplingParameter();
  const double lnG = log(Gamma);
  const double lnG2 = lnG*lnG;
  const double b0 = 0.258 - 0.0612*lnG + 0.0123*lnG2 - 1.0/Gamma;
  const double b1 = 0.0269 + 0.0318*lnG + 0.00814*lnG2;
  if (b0/b1 <= 0.0 || Gamma < 5.25 || Gamma > 171.8){
    const string msg = format<double>("Error: The IET schemes cannot be applied "
				      "to this state point because Gamma = %.8f "
				      "falls outside the range of validty of the "
				      "bridge function parameterization\n",
				      Gamma);
    throw runtime_error(msg); 
  }
  const double c1 = 0.498 - 0.280*lnG + 0.0294*lnG2;
  const double c2 = -0.412 + 0.219*lnG - 0.0251*lnG2;
  const double c3 = 0.0988 - 0.0534*lnG + 0.00682*lnG2;
  const double b02 = b0*b0;
  const double b03 = b02*b0;
  const double b04 = b03*b0;
  const double b05 = b04*b0;
  const double b06 = b05*b0;
  const double b07 = b06*b0;
  const double b08 = b07*b0;
  const double b12 = b1*b1;
  const double b13 = b12*b1;
  const double b14 = b13*b1;
  const double b15 = b14*b1;
  const double b16 = b15*b1;
  const double b17 = b16*b1;
  const double b18 = b17*b1;
  const double b02_b12 = b02/b12;
  const double b03_b13 = b03/b13;
  const double b04_b14 = b04/b14;
  const double b05_b15 = b05/b15;
  const double b06_b16 = b06/b16;
  const double b07_b17 = b07/b17; 
  const double b08_b18 = b08/b18;
  const double fact = sqrt(M_PI)/(4.0 * l2)*pow(b0/b1, 1.5);
  const double q2 = x*x;
  const double q3 = q2*x;
  const double q4 = q3*x;
  const double q5 = q4*x;
  const double q6 = q5*x;
  const double q7 = q6*x;
  const double q8 = q7*x;
  const double bf1 = -b0 + c1/16.0*(60.0*b02_b12 - 20.0*b03_b13*q2/l2
				    + b04_b14*q4/l4);
  const double bf2 = c2/64.0*(840.0*b03_b13 - 420.0*b04_b14*q2/l2 +
			      42.0*b05_b15*q4/l4 - b06_b16*q6/l6);;
  const double bf3 = c3/256.0*(15120.0*b04_b14 - 10080.0*b05_b15*q2/l2 +
			       1512.0*b06_b16*q4/l4 - 72.0*b07_b17*q6/l6 + 
			       b08_b18*q8/l8);
  return fact * q2 * (bf1 + bf2 + bf3) * exp(-b0*q2/(4.0*b1*l2));
}

double BridgeFunction::lct() const {
  const double Gamma = couplingParameter();
  auto func = [&](double r)->double{return lctIntegrand(r,Gamma);};
  itg.setR(x/lambda);
  itg.compute(func);
  return itg.getSolution() * (x/lambda) / Gamma;
  return 0.0;
}

double BridgeFunction::lctIntegrand(const double r, const double Gamma) const {
  if (Gamma < 5.0) {
    const string msg = format<double>("Error: The IET schemes cannot be applied "
				     "to this state point because Gamma = %.8f "
				     "falls outside the range of validty of the "
				     "bridge function parameterization\n",
				     Gamma);
   throw runtime_error(msg);
  }
  const double Gamma1_6 = pow(Gamma, 1./6.);
  const double lnG = log(Gamma);
  const double lnG2 = lnG*lnG; 
  const double lnG3 = lnG2*lnG;
  const double lnG4 = lnG3*lnG;
  const double a0 = Gamma * (0.076912 - 0.10465*lnG + 0.0056629*lnG2
  		       + 0.00025656*lnG3);
  const double a2 = Gamma * (0.068045 - 0.036952*lnG + 0.048818*lnG2
  		       - 0.0048985*lnG3);
  const double a3 = Gamma * (-0.30231 + 0.30457*lnG - 0.11424*lnG2
  		       + 0.0095993*lnG3);
  const double a4 = Gamma * (0.25111 - 0.26800*lnG + 0.082268*lnG2
  		       - 0.0064960*lnG3);
  const double a5 = Gamma * (-0.061894 + 0.066811*lnG - 0.019140*lnG2
  		       + 0.0014743*lnG3);
  const double c0 = Gamma * (0.25264 - 0.31615*lnG + 0.13135*lnG2 
		       - 0.023044*lnG3 + 0.0014666*lnG4);
  const double c1 = Gamma1_6 * (-12.665 + 20.802*lnG - 9.6296*lnG2 
			  + 1.7889*lnG3 - 0.11810*lnG4); 
  const double c2 = Gamma1_6 * (15.285 - 14.076*lnG + 5.7558*lnG2 
			  - 1.0188*lnG3 + 0.06551*lnG4); 
  const double c3 = Gamma1_6 * (35.330 - 40.727*lnG + 16.690*lnG2 
			  - 2.8905*lnG3 + 0.18243*lnG4);     
  const double r2 = r*r;
  const double r3 = r2*r;
  const double r4 = r3*r;
  const double r5 = r4*r;
  const double rshift = r - 1.44;
  const double bsr = a0 + a2*r2 + a3*r3 + a4*r4 + a5*r5;
  const double blr = c0 * exp(-c1*rshift) * exp(-0.3*r2)
    * ( cos(c2*rshift) + c3*exp(-3.5*rshift) );
  const double sf = 0.5 * ( 1.0 + erf(5.0*(r - 1.50)) );
  return r *( (1 - sf)*bsr + sf*blr);
}


