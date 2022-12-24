#include "stls.hpp"
#include "chemicalpotential.hpp"

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
  // Ideal density response
  cout << "Computing ideal density response: "; 
  computeIdr();
  cout << "Done." << endl;
  // Hartree-Fock static structure factor
  cout << "Computing HF static structure factor: "; 
  computeSsfHF();
  cout << "Done." << endl;
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
  Integrator1D itg;
  idr.resize(nl);
  for (int l=0; l<nl; ++l){
    computeIdrSingleFrequency(itg, l);
  }
}

void Stls::computeIdrSingleFrequency(Integrator1D &itg,
				     const int l){
  int nx = wvg.size();
  idr[l].resize(nx);
  for (int i=0; i<nx; ++i){
    Idr idrTmp(l, wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back());
    idr[l][i] = idrTmp.get(itg);
  }
}

// Compute Hartree-Fock static structure factor
void Stls::computeSsfHF(){
  assert(computedChemicalPotential);
  Integrator1D itg;
  int nx = wvg.size();
  ssfHF.resize(nx);
  for (int i=0; i<nx; ++i){
    SsfHF ssfTmp(wvg[i], in.getDegeneracy(), mu, wvg.front(), wvg.back());
    ssfHF[i] = ssfTmp.get(itg);
    cout << wvg[i] << " " << ssfHF[i] << endl;
  }
}


// -----------------------------------------------------------------
// Ideal density response class
// -----------------------------------------------------------------

// Integrand for ideal density response
double Idr::integrand(double y){
  if (l == 0)
    return x0(y);
  else
    return xl(y);
}

// Get at finite temperature
double Idr::get(Integrator1D &itg) {
  auto func = [this](double y)->double{return integrand(y);};
  itg.compute(func, yMin, yMax);
  return itg.getSolution();
}


// Integrand for frequency = l and wave-vector = x
double Idr::xl(double y) {
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
double Idr::x0(double y) {
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
double SsfHF::integrand(double y) {
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
double SsfHF::get(Integrator1D &itg) {
  auto func = [this](double y)->double{return integrand(y);};
  itg.compute(func, yMin, yMax);
  return 1.0 + itg.getSolution();
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
