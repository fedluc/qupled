#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <functional>
#include "util.hpp"

using namespace std;
using namespace inpututil;

#define NO_FILE_NAME ""

class StaticInput {

private:

  // Mixing parameter for the iterative procedure
  double aMix;
  // Minimum error for convergence in the iterative procedure
  double errMin;  
  // Wave-vector grid resolution
  double dx;
  // cutoff for the wave-vector grid
  double xmax;
  // Initial guess for the chemical potential calculation
  vector<double> muGuess;
  // Number of matsubara frequencies
  int nl;
  // Maximum number of iterations
  int nIter;
  // Output frequency
  int outIter;
  // Setters 
  void setMixingParameter(cString  &aMix);
  void setErrMin(cString &errMin);
  void setWaveVectorGridRes(cString &waveVectorGridRes);
  void setWaveVectorGridCutoff(cString  &waveVectorGridCutoff);
  void setChemicalPotentialGuess(cString &muGuessStr);
  void setNMatsubara(cString &nMatsubara);
  void setNIter(cString &nIter);
  void setOutIter(cString &outIter);
  // Helper methods to read the input file
  void assignInputToData(const string &keyword, const string &value);
  // Print content of the data structure
  void print() const;
  // Friends
  friend class Input;
  
public:
  
  //Constructor
  StaticInput();
  // Getters
  double getMixingParameter() const { return aMix; };
  double getErrMin() const { return errMin; };
  double getWaveVectorGridRes() const { return dx; };
  double getWaveVectorGridCutoff() const { return xmax; };
  vector<double> getChemicalPotentialGuess() const { return muGuess; };
  int getNMatsubara() const { return nl; };
  int getNIter() const { return nIter; };
  int getOutIter() const { return outIter; };
  
};

class StlsInput {

private:
 
  // Mapping between the quantum and classical state points for the IET-based schemes
  string IETMapping;
  // Name of the file used to store the restart data
  string restartFileName;
  // Setters
  void setIETMapping(cString &IETMapping);
  void setRestartFileName(cString &restartFileName);
  // Helper methods to read the input file
  void assignInputToData(const string &keyword, const string &value);
  // Print content of the data structure
  void print() const ;
  // Friends
  friend class Input;
  
public:

  //Constructor
  StlsInput();
  // Getters
  string getIETMapping() const { return IETMapping; };
  string getRestartFileName() const { return restartFileName; };


};

class QstlsInput { // Input properties for the qSTLS methods

private:

  // Use static approximation to compute the auxiliary density response (adr)
  bool useStaticAdr;
  // Name of the files used to store the fixed component of the adr
  string fixedFileName;
  // Setters 
  void setUseStaticAdr(cString &useStaticAdr);
  void setFixedFileName(cString &fixedFileName);
    // Helper methods to read the input file
  void assignInputToData(const string &keyword, const string &value);
   // Print content of the data structure
  void print() const ;
  // Friends
  friend class Input;

public:

  //Constructor
  QstlsInput();
  // Getters
  bool getUseStaticAdr() const { return useStaticAdr; };
  string getFixedFileName() const {return fixedFileName; };

};


class Input {

private:

  // theory to be solved
  string theory;
  // type of theory
  bool isClassicTheory;
  bool isQuantumTheory;
  // degeneracy parameter
  double Theta;
  // quantum coupling parameter
  double rs;
  // number of threads for parallel calculations
  size_t nThreads;
  // scheme for 2D integrals
  string int2DScheme;
  // input for static calcualtions
  StaticInput stat;
  // input for stls calculations
  StlsInput stls;
  // input for qstls calculations
  QstlsInput qstls;
  // Setters
  void setTheory(cString &theory);
  void setDegeneracy(cString &Theta);
  void setCoupling(cString &rs);
  void setThreads(cString &nThreads);
  void setInt2DScheme(cString &int2DScheme);
  // Helper methods to read the input file
  void parseInputLine(cString &line);
  void assignInputToData(cVector<string> &input);
  void assignInputToBaseData(cString &keyword, cString &value);
  void assignInputToStaticData(cString &keyword, cString &value);
  void assignInputToStlsData(cString &keyword, cString &value);
  void assignInputToQstlsData(cString &keyword, cString &value);
   
public:

  //Constructor
  Input();
  // Getters
  string getTheory() const { return theory; };
  bool isClassic() const { return isClassicTheory; };
  double getDegeneracy() const {return Theta; };
  double getCoupling() const { return rs; };
  size_t getNThreads() const { return nThreads; }
  string getInt2DScheme() const { return int2DScheme; }
  const StaticInput& getStaticInput() const { return stat; };
  const StlsInput& getStlsInput() const { return stls; };
  const QstlsInput& getQstlsInput() const { return qstls; };
  // Read input file
  void readInput(cString &fileName);
  void print() const;
  
};

// class vsInput { // Input properties for the vs-stls methods

// private:

//   // Solve compressibility sum-rule (csr)
//   bool solveCsr;
//   // Name of the file used to store the thermodynamic integration data
//   string thermoFileName;
//   // Initial guess for the free parameter
//   double alpha;
//   // Mixing parameter for the iterative procedure used to enforce the csr
//   double aMix;
//   // Resolution of the coupling parameter grid
//   double drs;
//   // Resolution of the degeneracy parameter grid
//   double dt;
//   // Minimum error for convergence in the iterations used to enforce the csr
//   double minErr;

// public:

//   //Constructor
//   vsInput();
//   // Getters
//   bool solveCsr();
//   string getThermoFileName();
//   double getAlpha();
//   double getAMix();
//   double getDrs();
//   double getDt();
//   double getMinErr();
//   // Setters
//   void SetSolveCsr(bool solveCsr);
//   void setThermoFileName(string thermoFileName);
//   void setAlpha(double alpha);
//   void setAMix(double aMix);
//   void setDrs(double drs);
//   void setDt(double dt);
//   void setMinErr(double minErr);

// };

// class dynInput { // Input properties for the dynamic properties

// private:

//   // Name of the file with the density response used for the calculation of the dynamic properties
//   string drFileName;
//   // Name of the file with the structural properties used for the calculation of the dynamic properties
//   string structFileName;
//   // Resolution of the frequency grid 
//   double dW;
//   // Lower cutoff of the frequency grid
//   double Wmin;
//   // Upper cutoff of the frequency grid
//   double Wmax;
//   // Wave-vector used for the output
//   double waveVector;

// public:

//   //Constructor
//   dynInput();
//   // Getters
//   string getDrFileName();
//   string getStructFileName();
//   double getDw();
//   double getWmin();
//   double getWmax();
//   double getWaveVector();
//   // Setters
//   void setDrFileName(string densityResponseFileName);
//   void setStructFileName(string structureFileName);
//   void setDw(double dW);
//   void setWmin(double Wmin);
//   void setWmax(double Wmax);
//   void setWaveVector(double waveVector);

// };

#endif
