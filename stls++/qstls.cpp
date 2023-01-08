#include <omp.h>
#include "util.hpp"
#include "qstls.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace binUtil;

// -----------------------------------------------------------------
// QSTLS class
// -----------------------------------------------------------------

void Qstls::compute(){
  if (in.getDegeneracy() == 0.0) {
    throw runtime_error("Ground state calculations are not available "
			"for the quantum schemes");
  }
  Stls::init();
  if (verbose) cout << "Structural properties calculation ..." << endl;
  doIterations();
  if (verbose) cout << "Done" << endl;
  if (verbose) cout << "Writing output files: ";
  writeOutput();
  writeRestart();
  if (verbose) cout << "Done" << endl;
}


// qstls iterations
void Qstls::doIterations() {
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int maxIter = statIn->getNIter();
  const int outIter = statIn->getOutIter();
  const double minErr = statIn->getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  while (counter < maxIter+1 && err > minErr ) {
    // Start timing
    double tic = omp_get_wtime();
    // Update auxiliary density response
    computeAdr();
    // Update static structure factor
    computeSsf();
    // Write output
    if (counter % outIter == 0) { writeOutput();};
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
}

// Compute auxiliary density response
void Qstls::computeAdr() {
  assert(itg != NULL);
  assert(ssfOld.size() > 0);
  const int nx = wvg.size();
  const int nl = in.getStaticInput()->getNMatsubara();
  if (adrFixed.size() == 0) computeAdrFixed();
  if (adr.size() == 0) adr.resize(nx);
  if (slfc.size() == 0) slfc.resize(nx);
  assert(wvg.size() == ssfOld.size());
  const shared_ptr<Interpolator> ssfi = make_shared<Interpolator>(wvg, ssfOld);
  for (int i=0; i<nx; ++i) {
    Adr adrTmp(nl, in.getDegeneracy(), wvg.front(),
	       wvg.back(), itg, ssfi);
    adrTmp.get(wvg, adrFixed[i], adr[i]);
    slfc[i] = adr[i][0];
  }
}

void Qstls::computeAdrFixed() {
  cout << "Computing fixed component of the auxiliary density response: ";
  fflush(stdout);
  assert(computedChemicalPotential);
  loadAdrFixed();
  if (adrFixed.size() > 0) { return; }
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int nx = wvg.size();
  const int nl = statIn->getNMatsubara();
  adrFixed.resize(nx);
  //#pragma omp for
  for (int i=0; i<nx; ++i) {
    const shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>();
    AdrFixed adrTmp(nl, wvg[i], in.getDegeneracy(), mu,
		    wvg.front(), wvg.back(), itg2);
    adrTmp.get(wvg, adrFixed[i]);
  }
  //writeRestart();
  cout << "Done" << endl;
}


void Qstls::loadAdrFixed() {
  const string fileName = in.getQstlsInput()->getFixedFileName();
  if (fileName == "") return;
  decltype(wvg) wvgFile;
  decltype(adrFixed) adrFixedFile;
  double ThetaFile;
  readRestart(fileName, wvgFile, adrFixedFile, ThetaFile);
  int nlFile = adrFixedFile[0].size();
  const double tol = 1e-15;
  const string errMsg = "Data loaded from file for the fixed component "
    "of the auxiliary density response is not "
    "compatible with input";
  if (nlFile != in.getStaticInput()->getNMatsubara()) {
    throw runtime_error(errMsg + ", wrong number of Matsubara frequencies.");
  }
  if (ThetaFile - in.getDegeneracy() > tol) {
    throw runtime_error(errMsg + ", wrong degeneracy parameter.");
  }
  if (wvgFile.size() != wvg.size()) {
    throw runtime_error(errMsg + ", wrong wave-vector grid.");
  }
  wvgFile = diff(wvgFile, wvg);
  if (*max_element(wvgFile.begin(), wvgFile.end()) > tol) {
    throw runtime_error(errMsg + ", wrong wave-vector grid.");
  }
  adrFixed = adrFixedFile;
}


// Compute static structure factor
void Qstls::computeSsf(){
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  computeSsfFinite();
}

// Compute static structure factor at finite temperature
void Qstls::computeSsfFinite(){
  assert(computedChemicalPotential);
  assert(adr.size() > 0);
  assert(idr.size() > 0);
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  if (ssf.size() == 0) ssf.resize(nx);
  for (int i=0; i<nx; ++i){
    Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], idr[i], adr[i]);
    ssf[i] = ssfTmp.get();
  }
}

// Initial guess for qstls iterations
void Qstls::initialGuess() {
  const int nx = wvg.size();
  ssfOld.resize(nx);
  ssf.resize(nx);
  // From file
  const string fileName = in.getStlsInput()->getRestartFileName();
  if (fileName.size() > 0) {
    vector<double> wvgFile;
    vector<double> ssfFile;
    readRestart(fileName, wvgFile, ssfFile);
    assert(wvgFile.size() == ssfFile.size());
    const Interpolator ssfi(wvgFile, ssfFile);
    const double xmaxi = wvgFile.back();
    for (int i=0; i<wvg.size(); ++i) {
      const double x = wvg[i];
      if (x <= xmaxi) { ssfOld[i] = ssfi.eval(x);}
      else { ssfOld[i] = 1.0; }
    }
    return;
  }
  // Default
  Stls stls(in, false, false);
  stls.compute();
  stls.getSsf(ssfOld);
}

// Compute residual error for the qstls iterations
double Qstls::computeError(){
  return rms(ssf, ssfOld, false);
}

// Update solution during qstls iterations
void Qstls::updateSolution(){
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const double aMix = statIn->getMixingParameter();
  ssfOld = sum(mult(ssf, aMix), mult(ssfOld, 1 - aMix));
}

// Write output files
void Qstls::writeOutput() const{
  Stls::writeOutput();
  writeAdr();
}

void Qstls::writeAdr() const {
  if (in.getDegeneracy() == 0.0) return;
  const string fileName = format<double,double>("adr_rs%.3f_theta%.3f_"
						+ in.getTheory() + ".dat",
						in.getCoupling(),
						in.getDegeneracy());
  ofstream file;
  file.open(fileName);
  if (!file.is_open()) {
    throw runtime_error("Output file " + fileName + " could not be created.");
  }
  for (int i=0; i<adr.size(); ++i){
    const string el1 = format<double>("%.8e", wvg[i]);
    file << el1;
    for (int l=0; l<adr[i].size(); ++l) {
      const string el2 = format<double>("%.8e ", adr[i][l]);
      file << el2;
    }
    file << endl;
  }
  file.close();
}

// Restart files
void Qstls::writeRestart() const {
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
  int nl = in.getStaticInput()->getNMatsubara();
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<int>(file, nl);
  writeDataToBinary<double>(file, in.getDegeneracy());
  writeDataToBinary<decltype(wvg)>(file, wvg);
  writeDataToBinary<decltype(ssf)>(file, ssf);
  writeDataToBinary<decltype(adrFixed)>(file, adrFixed);
  file.close();
  if (!file) {
    throw runtime_error("Error in writing the file " + fileName);
  }
}

void Qstls::readRestart(const string &fileName,
			decltype(wvg) &wvgFile,
			decltype(ssf) &ssfFile) const {
  decltype(adrFixed) tmp1;
  double tmp2;
  readRestart(fileName, wvgFile, ssfFile, tmp1, tmp2);
}

void Qstls::readRestart(const string &fileName,
			decltype(wvg) &wvgFile,
			decltype(adrFixed) &adrFixedFile,
			double &Theta) const {
  decltype(wvg) tmp;
  readRestart(fileName, wvgFile, tmp, adrFixedFile, Theta);
}

void Qstls::readRestart(const string &fileName,
			decltype(wvg) &wvgFile,
			decltype(ssf) &ssfFile,
			decltype(adrFixed) &adrFixedFile,
			double &Theta) const {
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Input file " + fileName + " could not be opened.");
  }
  int nx;
  int nl;
  readDataFromBinary<int>(file, nx);
  readDataFromBinary<int>(file, nl);
  readDataFromBinary<double>(file, Theta);
  if (!file) {
    throw runtime_error("Error in reading the file " + fileName);
  }
  cout << wvg.size() << " " << nx << " " << nl << " " << Theta << endl; 
  wvgFile.resize(nx);
  ssfFile.resize(nx);
  readDataFromBinary<decltype(wvgFile)>(file, wvgFile);
  readDataFromBinary<decltype(ssfFile)>(file, ssfFile);
  adrFixedFile.resize(nx);
  for (auto &el1 : adrFixedFile) {
    el1.resize(nl);
    for (auto &el2 : el1) {
      el2.resize(nx);
    }			 
  }
  readDataFromBinary<decltype(adrFixedFile)>(file, adrFixedFile);
  file.close();
}

// -----------------------------------------------------------------
// Auxiliary density response classes
// -----------------------------------------------------------------

// Compute static structure factor
double Adr::ssf(const double y) const {
  return ssfi->eval(y);
}
// Compute fixed component (from interpolator)
double Adr::fix(const double y) const {
  return fixi->eval(y);
}

// Integrand
double Adr::integrand(const double y) const{
  return y * fix(y) * (ssf(y) - 1.0);
}

// Get result of integration
void Adr::get(const vector<double> &wvg,
	      const vector<vector<double>> &fixed,
	      vector<double> &res) {
  res.resize(nl);
  for (int l = 0; l < nl; ++l){
    assert(fixed[l].size() > 0);
    assert(fixed[l].size() == wvg.size());
    fixi.reset(new Interpolator(wvg, fixed[l]));
    auto func = [&](double y)->double{return integrand(y);};
    itg->compute(func, yMin, yMax);
    res[l] = itg->getSolution();
    res[l] *= (l==0) ? -3.0/(4.0*Theta) : -3.0/8.0;
  }
}

// Get fixed component
void AdrFixed::get(vector<double> &wvg,
		   vector<vector<double>> &res) const {
  res.resize(nl);
  const int nx = wvg.size();
  for (int l = 0; l < nl; ++l){
    res[l].resize(nx);
    if (x == 0.0) {
      fill(res[l].begin(), res[l].end(), 0.0);
      continue;
    }
    for (int i = 0; i < nx; ++i) {
      auto tMin = [&](double q)->double{return x*x - x*wvg[i];};
      auto tMax = [&](double q)->double{return x*x + x*wvg[i];};
      auto func1 = [&](double q)->double{return integrand1(q, l);};
      auto func2 = [&](double t)->double{return integrand2(t, wvg[i], l);};
      itg->compute(func1, func2, qMin, qMax, tMin, tMax);
      res[l][i] = itg->getSolution();
    }
  }
}

double AdrFixed::integrand1(const double q, const double l) const {
  if (l == 0) return q/(exp(q*q/Theta - mu) + exp(-q*q/Theta + mu) + 2.0);
  return q/(exp(q*q/Theta - mu) + 1.0);
}

double AdrFixed::integrand2(const double t, const double y, const double l) const {
  const double q = itg->getX();
  if (q == 0) return 0;
  const double x2 = x*x;
  const double y2 = y*y;
  const double q2 = q*q;
  const double txq = 2.0 * x * q;
  if (l == 0) {
    if (t == txq) return 2.0*q2/(y2 + 2.0*txq - x2);
    const double t2 = t*t;
    double logarg = (t + txq)/(t - txq);
    logarg = (logarg < 0.0) ? -logarg : logarg;
    return 1.0/(2.0*t + y2 - x2)*((q2 - t2/(4.0*x2))*log(logarg) + q*t/x);
  }
  const double tplT = 2.0 * M_PI * l * Theta;
  const double tplT2 = tplT*tplT;
  const double txqpt = txq + t;
  const double txqmt = txq - t;
  const double txqpt2 = txqpt*txqpt;
  const double txqmt2 = txqmt*txqmt;
  const double logarg = (txqpt2 + tplT2)/(txqmt2 + tplT2);
  return 1.0/(2.0*t + y*y - x*x)*log(logarg);
}
