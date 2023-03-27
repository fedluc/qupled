#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include <iostream>

using namespace std;

#define EMPTY_STRING ""

class Input {

private:

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
  
};

class StlsInput : public Input {

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
  // Name of the file used to store the restart data
  string restartFileName;
  // cutoff for the wave-vector grid
  double xmax;
  
public:
  
  //Constructor
  StlsInput(const double &rs_,
	    const double &Theta_,
	    const string &theory_)
    : Input(rs_, Theta_, theory_), aMix(1.0), dx(0.1), errMin(1e-5),
      IETMapping("standard"), muGuess({-10, 10}), nl(128), nIter(1000),
      outIter(10), restartFileName(EMPTY_STRING), xmax(10.0) { ; };
  // Setters
  void setChemicalPotentialGuess(const double &muMin,
				 const double &muMax);
  void setErrMin(const double &errMin);
  void setMixingParameter(const double  &aMix);
  void setIETMapping(const string &IETMapping);
  void setNMatsubara(const int &nMatsubara);
  void setNIter(const int &nIter);
  void setOutIter(const int &outIter);
  void setRestartFileName(const string &restartFileName);
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
  string getRestartFileName() const { return restartFileName; };
  double getWaveVectorGridRes() const { return dx; };
  double getWaveVectorGridCutoff() const { return xmax; };
  // Print content of the data structure
  void print() const;
  
};


class QstlsInput {

private:

  // Name of the files used to store the fixed component of the auxiliary density response (adr)
  string fixedFileName;
  // Use static approximation to compute the adr
  bool useStaticAdr;

public:

  //Constructor
  QstlsInput()
    : fixedFileName(EMPTY_STRING), useStaticAdr(false) { ; };
  // Setters
  void setFixedFileName(const string &fixedFileName);
  void setUseStaticAdr(const bool &useStaticAdr);
  // Getters
  string getFixedFileName() const {return fixedFileName; };
  bool getUseStaticAdr() const { return useStaticAdr; };
  // Print content of the data structure
  void print() const ;

};

#endif
