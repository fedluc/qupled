#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <functional>

using namespace std;

#define NO_FILE_NAME "noFileName"

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
  size_t nl;
  // Maximum number of iterations
  size_t nIter;
  // Helper methods to read the input file
  // const string allowedKeywords[7] = {"mixing", "error", "waveVectorResolution",
  //                                    "waveVectorCutoff", "chemicalPotential",
  // 				     "matsubaraFrequencies", "iterations"};
  void assignInputToData(const string &keyword, const string &value);
  // Friends
  friend class Input;
  
public:
  
  //Constructor
  StaticInput();
  // Getters
  double getMixingParameter();
  double getErrMin();
  double getWaveVectorGridRes();
  double getWaveVectorGridCutoff();
  vector<double> getChemicalPotentialGuess();
  size_t getNMatsubara();
  size_t getNIter();
  // Setters 
  void setMixingParameter(const double aMix);
  void setErrMin(const double errMin);
  void setWaveVectorGridRes(const double waveVectorGridRes);
  void setWaveVectorGridCutoff(const double waveVectorGridCutoff);
  void setChemicalPotentialGuess(const string &muGuessStr);
  void setNMatsubara(const size_t nMatsubara);
  void setNIter(const size_t nIter);
  
};

class StlsInput {

private:
 
  // Mapping between the quantum and classical state points for the IET-based schemes
  string IETMapping;
  // Name of the file used to store the restart data
  string restartFileName;
  
public:

  //Constructor
  StlsInput();
  // Getters
  string getIETMapping();
  string getRestartFileName();
  // Setters 
  void setIETMapping(string IETMapping);
  void setRestartFileName(string restartFileName);

};


class Input {

private:

  // Types
  typedef const string cString;
  template<typename T> using cVector = const vector<T>;
  // theory to be solved
  string theory; 
  // degeneracy parameter
  double Theta;
  // quantum coupling parameter
  double rs;
  // Number of threads for parallel calculations
  size_t nThreads;
  // input for static calcualtions
  StaticInput stat;
  // input for stls calculations
  StlsInput stls;
  // Helper methods to read the input file
  static vector<string> tokenize(cString &str, const char separator);
  static void matchKeyAndData(cVector<string> &keyword,
			      cVector<string> &input,
			      map<string,function<void(cString&, cString&)>> &funcArr);
  void parseInputLine(cString &line);
  void assignInputToData(cVector<string> &input);
  void assignInputToBaseData(cString &keyword, cString &value);
  void assignInputToStaticData(cString &keyword, cString &value);
  void assignInputToStlsData(cString &keyword, cString &value);

  // Friends
  friend class StaticInput;
  friend class StlsInput;
  
  // --- Additional members that will be added later
  // string mode;
  // string ssfRestartFileName;
  // string slfcRestartFileName;
  // qstlsInput qstls;
  // vsInput vs;
  // dynInput dyn;
  
public:

  //Constructor
  Input();
  // Getters
  string getTheory();
  double getDegeneracy();
  double getCoupling();
  StaticInput getStatic();
  StlsInput getStsl();
  // Setters
  void setTheory(cString &theory);
  void setDegeneracy(cString &Theta);
  void setCoupling(cString &rs);
  void setThreads(cString &nThreads);
  void readInput(cString &fileName);
  
};
// class qstlsInput { // Input properties for the qSTLS methods

// private:

//   // Use static approximation to compute the auxiliary density response (adr)
//   bool useStaticAdr;
//   // Name of the file used to store the restart data
//   string restartFileName;
//   // Name of the files used to store the fixed component of the adr
//   string fixedFileName;
//   // Name of the files used to store the fixed component of the adr with qIET approximation
//   string fixedIETFileName;

// public:

//   //Constructor
//   qstlsInput();
//   // Getters
//   bool useStaticAdr();
//   string getRestartName();
//   string getFixedName(const bool iet);
//   // Setters 
//   void setUseStatic(bool useStaticAdr);
//   void setRestartName(string restartName, const bool iet);
//   void setFixedName(string fixedName);

// };

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
