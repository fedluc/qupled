#ifndef INPUT_HPP
#define INPUT_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include "util.hpp"

// -----------------------------------------------------------------
// Base class to handle input for the dielectric schemes
// -----------------------------------------------------------------

class Input {
    
protected:

  // Accuracy for the integrals
  double intError;
  // quantum coupling parameter
  double rs;
  // degeneracy parameter
  double Theta;
  // number of threads for parallel calculations
  int nThreads;
  // type of theory
  bool isClassicTheory;
  bool isQuantumTheory;
  // scheme for 2D integrals
  string int2DScheme;
  // theory to be solved
  string theory;

public:

  // Setters
  void setCoupling(const double &rs);
  void setDegeneracy(const double &Theta);
  void setInt2DScheme(const string &int2DScheme);
  void setIntError(const double &intError);
  void setNThreads(const int &nThreads);
  void setTheory(const string &theory);
  // Getters
  double getCoupling() const { return rs; }
  double getDegeneracy() const {return Theta; }
  string getInt2DScheme() const { return int2DScheme; }
  double getIntError() const { return intError; }
  int getNThreads() const { return nThreads; }
  string getTheory() const { return theory; }
  bool isClassic() const { return isClassicTheory; }
  // Print content of the data structure
  void print() const;
  // Compare two Input objects
  bool isEqual(const Input &in) const;
  
};

// -----------------------------------------------------------------
// Class to handle input for the random phase approximation
// -----------------------------------------------------------------

class RpaInput : public Input {
  
protected:

  // Wave-vector grid resolution
  double dx;
  // cutoff for the wave-vector grid
  double xmax;
  // Number of matsubara frequencies
  int nl;
  // Initial guess for the chemical potential calculation
  vector<double> muGuess;
  
public:
  
  // Setters
  void setChemicalPotentialGuess(const vector<double> &muGuess);
  void setNMatsubara(const int &nMatsubara);
  void setWaveVectorGridRes(const double &waveVectorGridRes);
  void setWaveVectorGridCutoff(const double  &waveVectorGridCutoff);
  // Getters
  vector<double> getChemicalPotentialGuess() const { return muGuess; }
  int getNMatsubara() const { return nl; }
  double getWaveVectorGridRes() const { return dx; }
  double getWaveVectorGridCutoff() const { return xmax; }
  // Print content of the data structure
  void print() const;
  // Compare two StlsInput objects
  bool isEqual(const RpaInput &in) const;
  
};

// -----------------------------------------------------------------
// Class to handle input for the STLS and STLS-IET schemes
// -----------------------------------------------------------------

class StlsInput : public RpaInput {

public:

  struct SlfcGuess {
    vector<double> wvg;
    vector<double> slfc;
    bool operator==(const SlfcGuess &other) const {
      return wvg == other.wvg && slfc == other.slfc;
    }
  };
  
protected:

  // Mixing parameter for the iterative procedure
  double aMix;
  // Minimum error for convergence in the iterative procedure
  double errMin;
  // Maximum number of iterations
  int nIter;
  // Output frequency
  int outIter;
  // Mapping between the quantum and classical state points for the IET-based schemes
  string IETMapping;
  // Name of the file used to store the recovery data
  string recoveryFileName;
  // Initial guess
  SlfcGuess guess;
  
public:
  
  // Setters
  void setErrMin(const double &errMin);
  void setMixingParameter(const double  &aMix);
  void setIETMapping(const string &IETMapping);
  void setNIter(const int &nIter);
  void setOutIter(const int &outIter);
  void setRecoveryFileName(const string &recoveryFileName);
  void setGuess(const SlfcGuess &guess);
  // Getters
  double getErrMin() const { return errMin; }
  string getIETMapping() const { return IETMapping; }
  double getMixingParameter() const { return aMix; }
  int getNIter() const { return nIter; }
  int getOutIter() const { return outIter; }
  string getRecoveryFileName() const { return recoveryFileName; }
  SlfcGuess getGuess() const { return guess; }
  // Print content of the data structure
  void print() const;
  // Compare two StlsInput objects
  bool isEqual(const StlsInput &in) const;
  
};

// -----------------------------------------------------------------
// Class to handle input for the VSSTLS
// -----------------------------------------------------------------

class VSStlsInput : public StlsInput {

public:

  struct FreeEnergyIntegrand {
    vector<double> grid;
    vector<vector<double>> integrand;
    bool operator==(const FreeEnergyIntegrand &other) const {
      return grid == other.grid && integrand == other.integrand;
    }
  };
  
private:

  // Initial guess for the free parameter
  vector<double> alphaGuess;
  // Resolution of the coupling parameter grid
  double drs;
  // Resolution of the degeneracy parameter grid
  double dTheta;
  // Minimum error for the iterations used to define the free parameter
  double errMinAlpha;
  // Maximum number of iterations used to define the free parameter
  int nIterAlpha;
  // Pre-computed free energy integrand
  FreeEnergyIntegrand fxcIntegrand;
  
public:

  // Setters
  void setAlphaGuess(const vector<double>  &alphaGuess);
  void setCouplingResolution(const double &drs);
  void setDegeneracyResolution(const double &dTheta);
  void setErrMinAlpha(const double &errMinAlpha);
  void setNIterAlpha(const int& nIterAlpha);
  void setFreeEnergyIntegrand(const FreeEnergyIntegrand &freeEnergyIntegrand);
  // Getters 
  vector<double> getAlphaGuess() const { return alphaGuess; }
  double getCouplingResolution() const { return drs; }
  double getDegeneracyResolution() const { return dTheta; }
  double getErrMinAlpha() const { return errMinAlpha; }
  double getNIterAlpha() const { return nIterAlpha; }
  FreeEnergyIntegrand getFreeEnergyIntegrand() const { return fxcIntegrand; }
  // Print content of the data structure
  void print() const;
  // Compare two VSStls objects
  bool isEqual( const VSStlsInput &in ) const;
  
};

// -----------------------------------------------------------------
// Class to handle input for the QSTLS and QSTLS-IET schemes
// -----------------------------------------------------------------

class QstlsInput : public StlsInput {

public:
  
  struct QstlsGuess {
    vector<double> wvg;
    vector<double> ssf;
    vecUtil::Vector2D adr;
    int matsubara = 0;
    bool operator==(const QstlsGuess &other) const {
      return wvg == other.wvg && ssf == other.ssf
	&& adr == other.adr && matsubara == other.matsubara;
    }
  };

private:

  // Name of the file with the fixed component of the auxiliary density response (adr)
  string fixed;
  // Name of the file with the fixed component of the adr for iet schemes
  string fixedIet;
  // Initial guess
  QstlsGuess guess;

public:

  // Setters
  void setFixed(const string &fixed);
  void setFixedIet(const string &fixedIet);
  void setGuess(const QstlsGuess &guess);
  // Getters
  string getFixed() const {return fixed; }
  string getFixedIet() const { return fixedIet; }
  QstlsGuess getGuess() const { return guess; }
  // Print content of the data structure
  void print() const ;
  // Compare two QstlsInput objects
   bool isEqual(const QstlsInput &in) const;
  
};

#endif
