#include <omp.h>
#include "util.hpp"
#include "qstls.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace binUtil;

// -----------------------------------------------------------------
// QSTLS class
// -----------------------------------------------------------------

int Qstls::compute(){
  try {
    // Set number of OMP threads
    omp_set_num_threads(in.getNThreads());
    // Solve scheme
    Stls::init();
    if (verbose) cout << "Structural properties calculation ..." << endl;
    doIterations();
    if (verbose) cout << "Done" << endl;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}


// qstls iterations
void Qstls::doIterations() {
  assert(false);
  // Throw error message for ground state calculations
  if (in.getDegeneracy() == 0.0) {
    throw runtime_error("Ground state calculations are not available "
			"for the quantum schemes");
  }
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
    // Update auxiliary density response
    computeAdr();
    // Update static structure factor
    computeSsf();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // Write output
    if (counter % outIter == 0) { writeRecovery(); };
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

// Initial guess for qstls iterations
void Qstls::initialGuess() {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  // Resize variables used in iterations
  ssf.resize(nx);
  adr.resize(nx, nl);
  ssfOld.resize(nx);
  if (useIet) { adrOld.resize(nx, nl); }
  // From recovery file
  const string fileName = in.getRecoveryFileName();
  if (fileName != EMPTY_STRING) {
    vector<double> wvg_;
    vector<double> ssf_;
    Vector2D adr_;
    Vector3D adrFixed_;
    double Theta;
    int nl_;
    readRecovery(fileName, wvg_, ssf_, adr_, adrFixed_, Theta, nl_);
    initialGuessSsf(wvg_, ssf_);
    if (useIet) { initialGuessAdr(wvg_, adr_); }
    if (checkAdrFixed(wvg_, Theta, nl_) == 0) {
      adrFixed = adrFixed_;
    }
    return;
  }
  // From guess in input
  if (qin.getGuess().wvg.size() > 0) {
    const auto &guess = qin.getGuess();
    initialGuessSsf(guess.wvg, guess.ssf);
    if (useIet) { initialGuessAdr(guess.wvg, guess.adr); }
    return;
  }  
  // Default
  Stls stls(in, false, false);
  int status = stls.compute();
  if (status != 0) {
    throw runtime_error("Failed to compute the default initial guess");
  }
  ssfOld = stls.getSsf();
  if (useIet) { adrOld.fill(0.0); }
}

void Qstls::initialGuessSsf(const vector<double> &wvg_,
			    const vector<double> &ssf_) {
  const int nx = wvg.size();
  const Interpolator1D ssfi(wvg_, ssf_);
  const double xMax = wvg_.back();
  for (int i=0; i<nx; ++i) {
    const double &x = wvg[i];
    ssfOld[i] = (x <= xMax) ? ssfi.eval(x) : 1.0;
  }
}

void Qstls::initialGuessAdr(const vector<double> &wvg_,
			    const Vector2D &adr_) {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const int nx_ = adr_.size(0);
  const int nl_ = adr_.size(1);
  if (nx_ == 0 || nl_ == 0) {
    adrOld.fill(0.0);
    return;
  }
  const double &xMax = wvg_.back();
  vector<Interpolator1D> itp(nl_);
  for (int l=0; l<nl_; ++l) {
    vector<double> tmp(nx_);
    for (int i = 0; i<nx_; ++i) {
      tmp[i] =  adr_(i,l);
    }
    itp[l].reset(wvg_[0], tmp[0], nx_);
  }
  for (int i=0; i<nx; ++i) {
    const double &x = wvg[i];
    if (x > xMax) {
      adrOld.fill(i, 0.0);
      continue;
    }
    for (int l=0; l<nl; ++l) {
      adrOld(i,l) = (l <= nl_) ? itp[l].eval(x) : 0.0;
    }
  }
}

// Compute auxiliary density response
void Qstls::computeAdr() {
  const int nx = wvg.size();
  assert(adr.size() > 0);
  assert(ssfOld.size() > 0);
  if (slfc.size() == 0) slfc.resize(nx);
  if (adrFixed.size() == 0) computeAdrFixed();
  const Interpolator1D ssfi(wvg, ssfOld);
#pragma omp parallel for
  for (int i=0; i<nx; ++i) {
    Integrator1D itgPrivate;
    Adr adrTmp(in.getDegeneracy(), wvg.front(),
	       wvg.back(), wvg[i], ssfi, itgPrivate);
    adrTmp.get(wvg, adrFixed, adr);
  }
  if (useIet) computeAdrIet();
  for (int i=0; i<nx; ++i) {slfc[i] = adr(i,0); };
}

// Compute static structure factor
void Qstls::computeSsf(){
  computeSsfFinite();
}

// Compute static structure factor at finite temperature
void Qstls::computeSsfFinite(){
  assert(computedChemicalPotential);
  assert(ssf.size() > 0);
  assert(adr.size() > 0);
  assert(idr.size() > 0);
  if (useIet) assert(bf.size() > 0);
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  const int nl = idr.size(1);
  for (int i=0; i<nx; ++i){
    const double bfi = (useIet) ? bf[i] : 0;
    Qssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], nl, &idr(i), &adr(i), bfi);
    ssf[i] = ssfTmp.get();
  }
}

// Compute residual error for the qstls iterations
double Qstls::computeError(){
  return rms(ssf, ssfOld, false);
}

// Update solution during qstls iterations
void Qstls::updateSolution(){
  const double aMix = in.getMixingParameter();
  ssfOld = sum(mult(ssf, aMix), mult(ssfOld, 1 - aMix));
  if (useIet) {
    Vector2D tmp = adr;
    adrOld.mult(1 - aMix);
    tmp.mult(aMix);
    adrOld.sum(tmp);
  }
}

void Qstls::computeAdrFixed() {
  // Check if it adrFixed can be loaded from input
  readAdrFixedFile(adrFixed, qin.getFixed(), false);
  if (adrFixed.size() > 0) { return; }
  // Compute from scratch
  cout << "Computing fixed component of the auxiliary density response: ";
  fflush(stdout);
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  adrFixed.resize(nx, nl, nx);
#pragma omp parallel for
  for (int i=0; i<nx; ++i) {
    Integrator2D itg2;
    AdrFixed adrTmp(in.getDegeneracy(), wvg.front(), wvg.back(),
		    wvg[i], mu, itgGrid, itg2);
    adrTmp.get(wvg, adrFixed);
  }
  // Write result to output file
  const string fileName = format<double,double>("adr_fixed_rs%.3f_theta%.3f_QSTLS.bin",
						in.getCoupling(), in.getDegeneracy());
  writeAdrFixedFile(adrFixed, fileName);
  cout << "Done" << endl;
}

void Qstls::writeAdrFixedFile(const Vector3D &res,
			      const string &fileName) const {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const double Theta = in.getDegeneracy();
  ofstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Output file " + fileName + " could not be created.");
  }
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<int>(file, nl);
  writeDataToBinary<double>(file, Theta);
  writeDataToBinary<vector<double>>(file, wvg);
  writeDataToBinary<Vector3D>(file, res);
  file.close();
  if (!file) {
    throw runtime_error("Error in writing to file " + fileName);
  }
}

void Qstls::readAdrFixedFile(Vector3D &res,
			     const string &fileName,
			     const bool iet) const {
  if (fileName == EMPTY_STRING) { return; }
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Input file " + fileName + " could not be opened.");
  }
  int nx_;
  int nl_;
  vector<double> wvg_;
  double Theta_;
  readDataFromBinary<int>(file, nx_);
  readDataFromBinary<int>(file, nl_);
  readDataFromBinary<double>(file, Theta_);
  wvg_.resize(nx);
  readDataFromBinary<vector<double>>(file, wvg_);
  if (iet) { res.resize(nl, nx, nx); }
  else { res.resize(nx, nl, nx); }
  readDataFromBinary<Vector3D>(file, res);
  file.close();
  if (!file) {
    throw runtime_error("Error in reading from file " + fileName);
  }
  if (checkAdrFixed(wvg_, Theta_, nl_) != 0) {
    throw runtime_error("Fixed component of the auxiliary density response"
			"loaded from from file is incompatible with input");
  }
}

int Qstls::checkAdrFixed(const vector<double> &wvg_,
				 const double Theta_,
				 const int nl_) const {
  constexpr double tol = 1e-15;
  const vector<double> wvgDiff = diff(wvg_, wvg);
  const double &wvgMaxDiff = abs(*max_element(wvgDiff.begin(), wvgDiff.end()));
  const bool consistentMatsubara = nl_ == in.getNMatsubara();
  const bool consistentTheta = abs(Theta_ - in.getDegeneracy()) <= tol;
  const bool consistentGrid = wvg_.size() == wvg.size() && wvgMaxDiff <= tol;
  if (!consistentMatsubara ||
      !consistentTheta ||
      !consistentGrid) {
    return 1;
  }
  return 0;
}

void Qstls::computeAdrIet() {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  assert(adrOld.size() > 0);
  // Compute bridge function
  if (bf.size() == 0) computeBf();
  // Compute fixed part
  computeAdrFixedIet();
  // Setup interpolators
  const Interpolator1D ssfi(wvg, ssfOld);
  const Interpolator1D bfi(wvg, bf);
  vector<Interpolator1D> dlfci(nl);
  Interpolator1D tmp(wvg, ssfOld);
  for (int l=0; l<nl; ++l) {
    vector<double> dlfc(nx);
    for (int i = 0; i<nx; ++i) {
      dlfc[i] = (idr(i,l) > 0.0) ? adrOld(i,l)/idr(i,l) : 0;
    }
    dlfci[l].reset(wvg[0], dlfc[0], nx);
  }
  // Compute qstls-iet contribution to the adr
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  Vector2D adrIet(nx, nl);
#pragma omp parallel for
  for (int i=0; i<nx; ++i) {
    Integrator2D itgPrivate;
    Vector3D adrFixedPrivate;
    readAdrFixedFile(adrFixedPrivate, adrFixedIetFileInfo.at(i).first, true);
    AdrIet adrTmp(in.getDegeneracy(), wvg.front(), wvg.back(),
		  wvg[i], ssfi, dlfci, bfi, itgGrid, itgPrivate);
    adrTmp.get(wvg, adrFixedPrivate, adrIet);
  }
  // Sum qstls and qstls-iet contributions to adr
  adr.sum(adrIet);
}

void Qstls::computeAdrFixedIet() {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  vector<int> idx;
  // Check which files have to be created
  getAdrFixedIetFileInfo();
  for (const auto &member : adrFixedIetFileInfo) {
    if(!member.second.second) { idx.push_back(member.first); };
  }
  if (idx.size() == 0) { return; }
  // Write necessary files
  cout << "Computing fixed component of the iet auxiliary density response: ";
  fflush(stdout);
#pragma omp parallel for
  for (int i=0; i<idx.size(); ++i) {
    Integrator1D itgPrivate;
    Vector3D res(nl, nx, nx);
    AdrFixedIet adrTmp(in.getDegeneracy(), wvg.front(), wvg.back(),
		       wvg[idx[i]], mu, itgPrivate);
    adrTmp.get(wvg, res);
    writeAdrFixedFile(res, adrFixedIetFileInfo.at(idx[i]).first);
  }
  // Update content of adrFixedIetFileInfo
  getAdrFixedIetFileInfo();
  cout << "Done" << endl;
}

void Qstls::getAdrFixedIetFileInfo() {
  const int nx = wvg.size();
  const double Theta = in.getDegeneracy();
  adrFixedIetFileInfo.clear();
  for (int i=0; i<nx; ++i) {
    string name = format<double,double>("adr_fixed_rs%.3f_theta%.3f_"
					+ in.getTheory() + "_wv%.5f.bin",
					in.getCoupling(), Theta, wvg[i]);
    if (qin.getFixedIet() != EMPTY_STRING) {
      __fs::filesystem::path fullPath = qin.getFixedIet();
      fullPath /= name;
      name = fullPath.string();
    }
    const bool found = __fs::filesystem::exists(name);
    const pair<string,bool> filePair = pair<string,bool>(name,found);
    adrFixedIetFileInfo.insert(pair<int,decltype(filePair)>(i,filePair));
  }
}

// Recovery files
void Qstls::writeRecovery() {
  ofstream file;
  file.open(recoveryFileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Recovery file " + recoveryFileName + " could not be created.");
  }
  int nx = wvg.size();
  int nl = in.getNMatsubara();
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<int>(file, nl);
  writeDataToBinary<double>(file, in.getDegeneracy());
  writeDataToBinary<vector<double>>(file, wvg);
  writeDataToBinary<vector<double>>(file, ssf);
  writeDataToBinary<Vector2D>(file, adr);
  writeDataToBinary<Vector3D>(file, adrFixed);
  file.close();
  if (!file) {
    throw runtime_error("Error in writing the recovery file " + recoveryFileName);
  }
}

void Qstls::readRecovery(const string &fileName,
			vector<double> &wvg_,
			Vector3D &adrFixed_,
			double &Theta,
			int &nl) const {
  vector<double> tmp1;
  Vector2D tmp2;
  readRecovery(fileName, wvg_, tmp1, tmp2, adrFixed_, Theta, nl);
}

void Qstls::readRecovery(const string &fileName,
			 vector<double> &wvg_,
			 vector<double> &ssf_,
			 Vector2D &adr_,
			 Vector3D &adrFixed_,
			 double &Theta,
			 int &nl) const {
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throw runtime_error("Input file " + fileName + " could not be opened.");
  }
  int nx;
  readDataFromBinary<int>(file, nx);
  readDataFromBinary<int>(file, nl);
  readDataFromBinary<double>(file, Theta);
  wvg_.resize(nx);
  ssf_.resize(nx);
  adr_.resize(nx, nl);
  adrFixed_.resize(nx, nl, nx);
  readDataFromBinary<vector<double>>(file, wvg_);
  readDataFromBinary<vector<double>>(file, ssf_);
  readDataFromBinary<Vector2D>(file, adr_);
  readDataFromBinary<Vector3D>(file, adrFixed_);
  file.close();
  if (!file) {
    throw runtime_error("Error in reading the file " + fileName);
  }
}

// -----------------------------------------------------------------
// Qssf class
// -----------------------------------------------------------------

// Get static structure factor
double Qssf::get() const {
  if (x == 0.0) return 0.0;
  const double f1 = 4.0*lambda*rs/M_PI;
  const double f2 = 1 - bf;
  const double x2 = x*x;
  double f3 = 0.0;
  for (int l=0; l<nl; ++l) {
    const double f4 = f2 * idr[l];
    const double f5 = 1.0 + f1/x2*(f4 - adr[l]);
    const double f6 = idr[l]*(f4 - adr[l])/f5;
    f3 += (l == 0) ? f6 : 2 * f6;
  }
  return ssfHF - 1.5 * f1/x2 * Theta * f3;
}


// -----------------------------------------------------------------
// AdrBase class
// -----------------------------------------------------------------

// Compute static structure factor
double AdrBase::ssf(const double y) const {
  return ssfi.eval(y);
}

// -----------------------------------------------------------------
// Adr class
// -----------------------------------------------------------------

// Compute fixed component
double Adr::fix(const double y) const {
  return fixi.eval(y);
}

// Integrand
double Adr::integrand(const double y) const{
  return y * fix(y) * (ssf(y) - 1.0);
}

// Get result of integration
void Adr::get(const vector<double> &wvg,
	      const Vector3D &fixed,
	      Vector2D &res) {
  const int nx = wvg.size();
  const int nl = fixed.size(1);
  auto it = lower_bound(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  if (x == 0.0) {
    res.fill(ix, 0.0);
    return;
  }
  for (int l = 0; l < nl; ++l){
    fixi.reset(wvg[0], fixed(ix,l), nx);
    auto func = [&](double y)->double{return integrand(y);};
    itg.compute(func, yMin, yMax);
    res(ix, l) = itg.getSolution();
    res(ix, l) *= (l==0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixed class
// -----------------------------------------------------------------

// Get fixed component
void AdrFixed::get(vector<double> &wvg,
		   Vector3D &res) const {
  const int nx = wvg.size();
  const int nl = res.size(1);
  if ( x == 0.0 ) { res.fill(0.0); };
  const double x2 = x*x;
  auto it = find(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  for (int l = 0; l < nl; ++l){
    for (int i = 0; i < nx; ++i) {
      const double xq = x*wvg[i];
      auto tMin = [&](double q)->double{return x2 - xq;};
      auto tMax = [&](double q)->double{return x2 + xq;};
      auto func1 = [&](double q)->double{return integrand1(q, l);};
      auto func2 = [&](double t)->double{return integrand2(t, wvg[i], l);};
      itg.compute(func1, func2, qMin, qMax, tMin, tMax, itgGrid);
      res(ix, l, i) = itg.getSolution();
    }
  }
}

// Integrands for the fixed component
double AdrFixed::integrand1(const double q,
			    const double l) const {
  if (l == 0) return q/(exp(q*q/Theta - mu) + exp(-q*q/Theta + mu) + 2.0);
  return q/(exp(q*q/Theta - mu) + 1.0);
}

double AdrFixed::integrand2(const double t,
			    const double y,
			    const double l) const {
  const double q = itg.getX();
  if (q == 0 || t == 0 || y == 0) { return 0; };
  const double x2 = x*x;
  const double y2 = y*y;
  const double q2 = q*q;
  const double txq = 2.0 * x * q;
  if (l == 0) {
    if (t == txq) { return 2.0*q2/(y2 + 2.0*txq - x2); };
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

// -----------------------------------------------------------------
// AdrIet class
// -----------------------------------------------------------------

// Compute dynamic local field correction
double AdrIet::dlfc(const double y,
		    const int l) const {
  return dlfci[l].eval(y);
}

// Compute auxiliary density response
double AdrIet::bf(const double y) const {
  return bfi.eval(y);
}

// Compute fixed component
double AdrIet::fix(const double x, const double y) const {
  return fixi.eval(x,y);
}

// Integrands
double AdrIet::integrand1(const double q,
			  const int l) const {
  if (q == 0.0) { return 0.0; }
  const double p1 = (1 - bf(q)) * ssf(q);
  const double p2 = dlfc(q,l) * (ssf(q) - 1.0);
  return (p1 -  p2 - 1.0) / q;
}

double AdrIet::integrand2(const double y) const {
  const double q = itg.getX();
  return y * fix(q,y) * (ssf(y) - 1.0);
}

// Get result of integration
void AdrIet::get(const vector<double> &wvg,
		 const Vector3D &fixed,
		 Vector2D &res) {
  const int nx = wvg.size();
  const int nl = fixed.size(0);
  auto it = lower_bound(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  if (x == 0.0) {
    res.fill(ix, 0.0);
    return;
  }
  for (int l = 0; l < nl; ++l) {
    fixi.reset(wvg[0], wvg[0], fixed(l), nx, nx);
    auto yMin = [&](double q)->double{return (q > x) ? q - x : x - q;};
    auto yMax = [&](double q)->double{return min(qMax, q + x);};
    auto func1 = [&](double q)->double{return integrand1(q, l);};
    auto func2 = [&](double y)->double{return integrand2(y);};
    itg.compute(func1, func2, qMin, qMax, yMin, yMax, itgGrid);
    res(ix, l) = itg.getSolution();
    res(ix, l) *= (l == 0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixedIet class
// -----------------------------------------------------------------

// get fixed component
void AdrFixedIet::get(vector<double> &wvg,
		      Vector3D &res) const {
  if (x == 0) {
    res.fill(0.0);
    return;
  }
  const int nx = wvg.size();
  const int nl = res.size(0);
  for (int l = 0; l < nl; ++l){
    for (int i = 0; i < nx; ++i) {
      if (wvg[i] == 0.0) {
	res.fill(l, i, 0.0);
	continue;
      }
      for (int j = 0; j < nx; ++j) {
	auto func = [&](double t)->double{
	  return integrand(t, wvg[j], wvg[i], l);
	};
	itg.compute(func, tMin, tMax);
	res(l,i,j) = itg.getSolution();
      }  
    }
  }
}

// Integrand for the fixed component
double AdrFixedIet::integrand(const double t, const double y,
			      const double q, const double l) const {
  // l ---> l
  // x ---> x
  // u ---> q
  // w ---> y
  // y ---> t
  const double x2 = x*x;
  const double q2 = q*q;
  const double y2 = y*y;
  const double t2 = t*t;
  const double fxt = 4.0*x*t;
  const double qmypx = q2 - y2 + x2;
  if (l == 0) {
    double logarg = (qmypx + fxt)/(qmypx - fxt);
    if (logarg < 0.0) logarg = -logarg;
    return t / (exp(t2/Theta - mu) + exp(-t2/Theta + mu) + 2.0)*
      ((t2 - qmypx*qmypx/(16.0*x2))*log(logarg) + (t/x)*qmypx/2.0);
  }
  const double fplT = 4.0*M_PI*l*Theta;
  const double fplT2 = fplT*fplT;
  const double logarg = ((qmypx + fxt)*(qmypx + fxt) + fplT2)/
    ((qmypx - fxt)*(qmypx - fxt) + fplT2);
  return t / (exp(t2/Theta - mu) + 1.0)*log(logarg);
}