#ifndef QSTLS_HPP
#define QSTLS_HPP

#include <vector>
#include "input.hpp"
#include "stls.hpp"
#include "numerics.hpp"

using namespace std;
using namespace vecUtil;

// -----------------------------------------------------------------
// Solver for the STLS-based schemes
// -----------------------------------------------------------------

class Qstls : public Stls {
  
private: 

  // Auxiliary density response
  Vector2D<double> adr;
  Vector2D<double> adrOld;
  Vector3D<double> adrFixed;
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

// -----------------------------------------------------------------
// Classes for the auxiliary density response
// -----------------------------------------------------------------

class AdrBase {

protected:
  
  // Number of matsubara frequency
  const int nl;
  // Degeneracy parameter
  const double Theta;
  // Integration limits
  const double yMin;
  const double yMax;
  // Interpolator for the static structure factor
  const Interpolator &ssfi;
  // Interpolator for the fixed component
  Interpolator fixi;
  // Compute static structure factor
  double ssf(const double y) const;
  
public:

  // Constructor
  AdrBase(const int nl_,
	  const double Theta_,
	  const double yMin_,
	  const double yMax_,
	  const Interpolator &ssfi_) : nl(nl_), Theta(Theta_), yMin(yMin_),
				       yMax(yMax_), ssfi(ssfi_) {;};
  
};

class AdrFixedBase {

protected:

  // Number of matsubara frequency
  const int nl;
  // Degeneracy parameter
  const double Theta;
  // Integration limits
  const double qMin;
  const double qMax;
  // Wave-vector
  const double x;
  // Chemical potential
  const double mu;
  
public:

  // Constructor for finite temperature calculations
  AdrFixedBase(const int nl_,
	       const double Theta_,
	       const double qMin_,
	       const double qMax_,
	       const double x_,
	       const double mu_) : nl(nl_), Theta(Theta_), qMin(qMin_),
				   qMax(qMax_), x(x_), mu(mu_) {;};
  
};

class Adr : public AdrBase {

private:

  // Compute fixed component
  double fix(const double y) const;
  // integrand 
  double integrand(const double y) const;
  // Integrator object
  Integrator1D &itg;
  
public:

  // Constructor for finite temperature calculations
  Adr(const int nl_,
      const double Theta_,
      const double yMin_,
      const double yMax_,
      const Interpolator &ssfi_,
      Integrator1D &itg_) : AdrBase(nl_, Theta_, yMin_, yMax_, ssfi_),
			    itg(itg_) {;};
  // Get result of integration
  void get(const vector<double> &wvg,
	   const Vector2D<double> &fixed,
	   vector<double> &res);
  
};

class AdrFixed : public AdrFixedBase {
  
private:
  
  // Integrands 
  double integrand1(const double q,
		    const double l) const;
  double integrand2(const double t,
		    const double y,
		    const double l) const;
  // Integrator object
  Integrator2D &itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixed(const int nl_,
	   const double Theta_,
	   const double qMin_,
	   const double qMax_,
	   const double x_,
	   const double mu_,
	   Integrator2D &itg_) : AdrFixedBase(nl_, Theta_, qMin_, qMax_, x_, mu_),
				 itg(itg_) {;};
  // Get integration result
  void get(vector<double> &wvg,
	   Vector2D<double> &res) const;
  
};

// Class for the auxiliary density response calculation in the IET scheme
class AdrIet : public AdrBase {

private:

  // Wave-vector
  const double x;
  // Matsubara frequency
  const int l;
   // Integration limits
  const double &qMin = yMin;
  const double &qMax = yMax;
  // Integrands 
  double integrand1(const double q) const;
  double integrand2(const double y) const;
  // Integrator object
  Integrator2D &itg;
  // Interpolator for the ideal density response
  const Interpolator &idri;
  // Interpolator for the auxiliary density response
  const Interpolator &adri;
  // Interpolator for the bridge function contribution
  const Interpolator &bfi;
  // Interpolator for the fixed component 
  Interpolator2D fixi;
  // Compute ideal density response
  double idr(const double y) const;
  // Compute auxiliary density
  double adr(const double y) const;
  // Compute bridge function contribution
  double bf(const double y) const;
  // Compute fixed component
  double fix(const double x, const double y) const;
  
public:

  // Constructor for finite temperature calculations
  AdrIet(const double Theta_,
	 const double qMin_,
	 const double qMax_,
	 const double x_,
	 const int l_,
	 const Interpolator &ssfi_,
	 const Interpolator &idri_,
	 const Interpolator &adri_,
	 const Interpolator &bfi_,
	 Integrator2D &itg_) : AdrBase(0, Theta_, qMin_, qMax_, ssfi_),
			       x(x_), l(l_), itg(itg_), idri(idri_),
			       adri(adri_), bfi(bfi_) {;};
  // Get integration result
  void get(const vector<double> &wvg,
	   const Vector2D<double> &fixed,
	   double &res);
  
};

class AdrFixedIet : public AdrFixedBase {

private:

  // Integration limits
  const double &tMin = qMin;
  const double &tMax = qMax;
  // Integrands 
  double integrand(const double t, const double y,
		   const double q, const double l) const;
  // Integrator object
  Integrator1D &itg;
  
public:

  // Constructor for finite temperature calculations
  AdrFixedIet(const int nl_,
	      const double Theta_,
	      const double qMin_,
	      const double qMax_,
	      const double x_,
	      const double mu_,
	      Integrator1D &itg_) : AdrFixedBase(nl_, Theta_, qMin_, qMax_, x_, mu_),
				    itg(itg_) {;};
  // Get integration result
  void get(vector<double> &wvg,
	   Vector3D<double> &res) const;
  
};


#endif
