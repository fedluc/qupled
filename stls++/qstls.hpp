#ifndef QSTLS_HPP
#define QSTLS_HPP

#include <vector>
#include "input.hpp"
#include "stls.hpp"
#include "numerics.hpp"

using namespace std;

class Qstls : public Stls {
  
private: 

  // Auxiliary density response
  vector<vector<double>> adr;
  vector<vector<vector<double>>> adrFixed;
  // Static structure factor (for iterations)
  vector<double> ssfOld;
  // Compute auxiliary density response
  void computeAdr();
  void computeAdrFixed();
  void loadAdrFixed();
  // Compute static structure factor at finite temperature
  void computeSsf();
  void computeSsfFinite();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  double computeError();
  void updateSolution();
   // Write output files
  void writeOutput() const;
  void writeAdr() const;
  // Restart files
  void writeRestart() const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvgFile,
		   decltype(ssf) &ssfFile) const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvgFile,
		   decltype(adrFixed) &adrFixedFile,
		   double &Theta) const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvgFile,
		   decltype(ssf) &ssfFile,
		   decltype(adrFixed) &adrFixedFile,
		   double &Theta) const;

public:

  // Constructors
  Qstls(const Input in_) : Stls(in_) {;};
  // Compute qstls scheme
  void compute(); 

};

// Class for the auxiliary density response calculation
class Adr {

private:
  
  // Number of matsubara frequency
  const int nl = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Integration limits
  const double yMin = 0;
  const double yMax = 0;
  // integrand 
  double integrand(const double y) const;
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  // Interpolator for the static structure factor
  const shared_ptr<Interpolator> ssfi;
  // Interpolator for the fixed component
  shared_ptr<Interpolator> fixi;
  // Compute static structure factor
  double ssf(const double y) const;
  // Compute fixed component
  double fix(const double y) const;
  
public:

  // Constructor for finite temperature calculations
  Adr(const int nl_,
      const double Theta_,
      const double yMin_,
      const double yMax_,
      const shared_ptr<Integrator1D> &itg_,
      const shared_ptr<Interpolator> &ssfi_)
    : nl(nl_), Theta(Theta_), yMin(yMin_),
      yMax(yMax_), itg(itg_), ssfi(ssfi_) {;};
  // Get result of integration
  void get(const vector<double> &wvg,
	   const vector<vector<double>> &fixed,
	   vector<double> &res);
  
};

class AdrFixed {

private:

  // Number of matsubara frequencies
  const int nl = 0;
  // Wave-vectors
  const double x = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Chemical potential
  const double mu = 0;
  // Integration limits
  const double qMin = 0;
  const double qMax = 0;
  // Integrands 
  double integrand1(const double q, const double l) const;
  double integrand2(const double t, const double y, const double l) const;
  // Integrator object
  const shared_ptr<Integrator2D> itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixed(const int nl_,
	   const double x_,
	   const double Theta_,
	   const double mu_,
	   const double qMin_,
	   const double qMax_,
	   const shared_ptr<Integrator2D> &itg_)
    : nl(nl_), x(x_), Theta(Theta_), mu(mu_),
      qMin(qMin_), qMax(qMax_), itg(itg_) {;};
  // Get integration result
  void get(vector<double> &wvg,
	   vector<vector<double>> &res) const;
  
};


class AdrFixedIet {

private:

  // Number of matsubara frequencies
  const int nl = 0;
  // Wave-vectors
  const double x = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Chemical potential
  const double mu = 0;
  // Integration limits
  const double tMin = 0;
  const double tMax = 0;
  // Integrands 
  double integrand(const double t, const double y,
		   const double q, const double l) const;
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixedIet(const int nl_,
	      const double x_,
	      const double Theta_,
	      const double mu_,
	      const double tMin_,
	      const double tMax_,
	      const shared_ptr<Integrator1D> &itg_)
    : nl(nl_), x(x_), Theta(Theta_), mu(mu_),
      tMin(tMin_), tMax(tMax_), itg(itg_) {;};
  // Get integration result
  void get(vector<double> &wvg,
	   vector<vector<vector<double>>> &res) const;
  
};


#endif
