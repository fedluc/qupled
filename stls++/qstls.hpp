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
  vector<vector<double>> adrOld;
  vector<vector<vector<double>>> adrFixed;
  map<int,pair<string,bool>> adrFixedIetFileInfo;
  // Static structure factor (for iterations)
  vector<double> ssfOld;
  // Compute auxiliary density response
  void computeAdr();
  void computeAdrFixed();
  void loadAdrFixed();
  void checkAdrFixedFromFile(const decltype(wvg) &wvg_,
			     const double Theta_,
			     const int nl_) const;
  void computeAdrIet();
  void computeAdrFixedIet();
  void getAdrFixedIetFileInfo();
  void writeAdrFixedIetFile(const decltype(adrFixed) &res,
			    const int i) const;
  void readAdrFixedIetFile(decltype(adrFixed) &res,
			   const int i) const;
  // Compute static structure factor at finite temperature
  void computeSsf();
  void computeSsfFinite();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  void initialGuessSsf(const decltype(wvg) &wvg_,
		       const decltype(ssf) &adr_);
  void initialGuessAdr(const decltype(wvg) &wvg_,
		       const decltype(adr) &adr_);
  double computeError();
  void updateSolution();
   // Write output files
  void writeOutput() const;
  void writeAdr() const;
  // Restart files
  void writeRestart() const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvg_,
		   decltype(ssf) &ssf_,
		   decltype(adr) &adr_) const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvg_,
		   decltype(adrFixed) &adrFixed_,
		   double &Theta) const;
  void readRestart(const string &fileName,
		   decltype(wvg) &wvg_,
		   decltype(ssf) &ssf_,
		   decltype(adr) &adr,
		   decltype(adrFixed) &adrFixed_,
		   double &Theta) const;
  // Check if iet schemes should be used
  void checkIet() { useIet = in.getTheory() == "QSTLS-HNC" ||
      in.getTheory() == "QSTLS-IOI" ||
      in.getTheory() == "QSTLS-LCT";}

public:

  // Constructors
  Qstls(const Input in_) : Stls(in_) {checkIet();};
  // Compute qstls scheme
  void compute(); 

};

// Class for the auxiliary density response calculation
class Adr {

protected:
  
  // Number of matsubara frequency
  const int nl;
  // Degeneracy parameter
  const double Theta;
  // Integration limits
  const double yMin;
  const double yMax;

private:
  
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

class AdrFixed : public Adr {

protected:

  // Wave-vector
  const double x;
  // Chemical potential
  const double mu;

private:
  
  // Integration limits
  const double &qMin = yMin;
  const double &qMax = yMax;
  // Integrands 
  double integrand1(const double q, const double l) const;
  double integrand2(const double t, const double y, const double l) const;
  // Integrator object
  const shared_ptr<Integrator2D> itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixed(const int nl_,
	   const double Theta_,
	   const double qMin_,
	   const double qMax_,
	   const double x_,
	   const double mu_,
	   const shared_ptr<Integrator2D> &itg_)
    : Adr(nl_, Theta_, qMin_, qMax_, NULL, NULL),
      x(x_), mu(mu_), itg(itg_) {;};
  // Get integration result
  void get(vector<double> &wvg,
	   vector<vector<double>> &res) const;
  
};


class AdrFixedIet : public AdrFixed {

private:

  // Integration limits
  const double &tMin = yMin;
  const double &tMax = yMax;
  // Integrands 
  double integrand(const double t, const double y,
		   const double q, const double l) const;
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixedIet(const int nl_,
	      const double Theta_,
	      const double x_,
	      const double mu_,
	      const double tMin_,
	      const double tMax_,
	      const shared_ptr<Integrator1D> &itg_)
    : AdrFixed(nl_, Theta_, x_, mu_, tMin_, tMax_, NULL), itg(itg_){;};
  // Get integration result
  void get(vector<double> &wvg,
	   vector<vector<vector<double>>> &res) const;
  
};


#endif
