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
  if (!computedChemicalPotential) {
    computeChemicalPotential();
  }
  const shared_ptr<StaticInput> &statIn = in.getStaticInput();
  const int nl = statIn->getNMatsubara();
  Integrator1D itg;
  // Loop over the Mastubara frequencies
  idr.resize(nl);
  for (int l=0; l<nl; ++l){
    computeIdrSingleFrequency(itg, l);
  }
}

void Stls::computeIdrSingleFrequency(Integrator1D &itg,
				     const int l){
  int nx = wvg.size();
  idr[l].resize(nx);
  // Loop over all wave-vector grid elements
  for (int i=0; i<nx; ++i){
    Idr idri(l, wvg[i], in.getDegeneracy(), mu);
    auto func = [&idri](double y)->double{return idri.get(y);};
    itg.compute(func, wvg.front(), wvg.back());
    idr[l][i] = itg.getSolution();
    cout << idr[l][i] << endl;
  }
}

// -----------------------------------------------------------------
// Ideal density response class
// -----------------------------------------------------------------

// Get integrand for ideal density response
double Idr::get(double y){
  if (l == 0)
    return x0(y);
  else
    return xl(y);
}

// Partial ideal density response (frequency = l, vector = x)
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

// Partial ideal density response (frequency = 0, vector = x)
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

// -----------------------------------------------------------------
// Static structure factor class
// -----------------------------------------------------------------

