#ifndef INPUT_HPP
#define INPUT_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include "util.hpp"

#define EMPTY_STRING ""

class Input {

protected:

  // scheme for 2D integrals
  string int2DScheme;
  // type of theory
  bool isClassicTheory;
  bool isQuantumTheory;
  // number of threads for parallel calculations
  int nThreads;
  // quantum coupling parameter
  double rs;
  // theory to be solved
  string theory;
  // degeneracy parameter
  double Theta;
  // Initializers
  double initCoupling(const double &rs_);
  double initDegeneracy(const double &Theta_);
  string initTheory(const string &theory_);
  
public:

  //Constructor
  Input(const double &rs_,
	const double &Theta_,
	const string &theory_)
    : int2DScheme("full"), isClassicTheory(false), isQuantumTheory(false),
      nThreads(1), rs(initCoupling(rs_)), theory(initTheory(theory_)),
      Theta(initDegeneracy(Theta_)) { ; };
  // Setters
  void setCoupling(const double &rs);
  void setDegeneracy(const double &Theta);
  void setInt2DScheme(const string &int2DScheme);
  void setNThreads(const int &nThreads);
  void setTheory(const string &theory);
  // Getters
  double getCoupling() const { return rs; };
  double getDegeneracy() const {return Theta; };
  string getInt2DScheme() const { return int2DScheme; }
  int getNThreads() const { return nThreads; }
  string getTheory() const { return theory; };
  bool isClassic() const { return isClassicTheory; };
  // Print content of the data structure
  void print() const;
  // Compare two Input objects
  bool isEqual(const Input &in) const;
  
};

class StlsInput : public Input {

public:
  
  struct SlfcGuess {
    vector<double> wvg = vector<double>(0);
    vector<double> slfc = vector<double>(0);
  };
  
private:

  // Mixing parameter for the iterative procedure
  double aMix;
  // Wave-vector grid resolution
  double dx;
  // Minimum error for convergence in the iterative procedure
  double errMin;
  // Mapping between the quantum and classical state points for the IET-based schemes
  string IETMapping;
  // Initial guess for the chemical potential calculation
  vector<double> muGuess;
  // Number of matsubara frequencies
  int nl;
  // Maximum number of iterations
  int nIter;
  // Output frequency
  int outIter;
  // Name of the file used to store the recovery data
  string recoveryFileName;
  // Initial guess
  SlfcGuess guess;
  // cutoff for the wave-vector grid
  double xmax;
  
public:
  
  //Constructor
  StlsInput(const double &rs_,
	    const double &Theta_,
	    const string &theory_)
    : Input(rs_, Theta_, theory_), aMix(1.0), dx(0.1), errMin(1e-5),
      IETMapping("standard"), muGuess({-10, 10}), nl(128), nIter(1000),
      outIter(10), recoveryFileName(EMPTY_STRING), xmax(10.0) { ; };
  // Setters
  void setChemicalPotentialGuess(const vector<double> &muGuess);
  void setErrMin(const double &errMin);
  void setMixingParameter(const double  &aMix);
  void setIETMapping(const string &IETMapping);
  void setNMatsubara(const int &nMatsubara);
  void setNIter(const int &nIter);
  void setOutIter(const int &outIter);
  void setRecoveryFileName(const string &recoveryFileName);
  void setGuess(const SlfcGuess &guess);
  void setWaveVectorGridRes(const double &waveVectorGridRes);
  void setWaveVectorGridCutoff(const double  &waveVectorGridCutoff);
  // Getters
  vector<double> getChemicalPotentialGuess() const { return muGuess; };
  double getErrMin() const { return errMin; };
  string getIETMapping() const { return IETMapping; };
  double getMixingParameter() const { return aMix; };
  int getNMatsubara() const { return nl; };
  int getNIter() const { return nIter; };
  int getOutIter() const { return outIter; };
  string getRecoveryFileName() const { return recoveryFileName; };
  SlfcGuess getGuess() const { return guess; };
  double getWaveVectorGridRes() const { return dx; };
  double getWaveVectorGridCutoff() const { return xmax; };
  // Print content of the data structure
  void print() const;
  // Compare two StlsInput objects
  bool isEqual(const StlsInput &in) const;
  
};


class QstlsInput {

public:
  
  struct QstlsGuess {
    vector<double> wvg = vector<double>(0);
    vector<double> ssf = vector<double>(0);
    vecUtil::Vector2D adr;
    int matsubara = 0;
  };

private:

  // Name of the files used to store the fixed component of the auxiliary density response (adr)
  string fixed;
  string fixedIet;
  // Initial guess
  QstlsGuess guess;

public:

  //Constructor
  QstlsInput()
    : fixed(EMPTY_STRING), fixedIet(EMPTY_STRING) { ; };
  // Setters
  void setFixed(const string &fixed);
  void setFixedIet(const string &fixedIet);
  void setGuess(const QstlsGuess &guess);
  // Getters
  string getFixed() const {return fixed; };
  string getFixedIet() const { return fixedIet; };
  QstlsGuess getGuess() const { return guess; };
  // Print content of the data structure
  void print() const ;
  // Compare two QstlsInput objects
  bool isEqual(const QstlsInput &in) const;
  
};

#endif
