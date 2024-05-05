#ifndef INPUT_HPP
#define INPUT_HPP

#include <cassert>
#include <iostream>
#include <vector>

// Forward declarations
namespace vecUtil {
  class Vector2D;
}

// -----------------------------------------------------------------
// Base class to handle input for the dielectric schemes
// -----------------------------------------------------------------

class Input {

public:

  // Constructor
  explicit Input()
      : intError(0),
        rs(0),
        Theta(0),
        nThreads(0),
        isClassicTheory(false),
        isQuantumTheory(false),
        int2DScheme(""),
        theory("") {}
  // Setters
  void setCoupling(const double &rs);
  void setDegeneracy(const double &Theta);
  void setInt2DScheme(const std::string &int2DScheme);
  void setIntError(const double &intError);
  void setNThreads(const int &nThreads);
  void setTheory(const std::string &theory);
  // Getters
  double getCoupling() const { return rs; }
  double getDegeneracy() const { return Theta; }
  std::string getInt2DScheme() const { return int2DScheme; }
  double getIntError() const { return intError; }
  int getNThreads() const { return nThreads; }
  std::string getTheory() const { return theory; }
  bool isClassic() const { return isClassicTheory; }
  // Print content of the data structure
  void print() const;
  // Compare two Input objects
  bool isEqual(const Input &in) const;

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
  std::string int2DScheme;
  // theory to be solved
  std::string theory;
};

// -----------------------------------------------------------------
// Class to handle input for the random phase approximation
// -----------------------------------------------------------------

class RpaInput : public Input {

public:

  // Constructor
  explicit RpaInput()
      : dx(0),
        xmax(0),
        nl(0),
        muGuess(std::vector<double>(2, 0)) {}
  // Setters
  void setChemicalPotentialGuess(const std::vector<double> &muGuess);
  void setNMatsubara(const int &nMatsubara);
  void setWaveVectorGridRes(const double &waveVectorGridRes);
  void setWaveVectorGridCutoff(const double &waveVectorGridCutoff);
  // Getters
  std::vector<double> getChemicalPotentialGuess() const { return muGuess; }
  int getNMatsubara() const { return nl; }
  double getWaveVectorGridRes() const { return dx; }
  double getWaveVectorGridCutoff() const { return xmax; }
  // Print content of the data structure
  void print() const;
  // Compare two StlsInput objects
  bool isEqual(const RpaInput &in) const;

protected:

  // Wave-vector grid resolution
  double dx;
  // cutoff for the wave-vector grid
  double xmax;
  // Number of matsubara frequencies
  int nl;
  // Initial guess for the chemical potential calculation
  std::vector<double> muGuess;
};

// -----------------------------------------------------------------
// Class to handle input for the STLS and STLS-IET schemes
// -----------------------------------------------------------------

class StlsInput : public RpaInput {

public:

  // Typedef
  struct SlfcGuess {
    std::vector<double> wvg;
    std::vector<double> slfc;
    bool operator==(const SlfcGuess &other) const {
      return wvg == other.wvg && slfc == other.slfc;
    }
  };
  // Contructor
  explicit StlsInput()
      : aMix(0),
        errMin(0),
        nIter(0),
        outIter(0),
        IETMapping(""),
        recoveryFileName("") {}
  // Setters
  void setErrMin(const double &errMin);
  void setMixingParameter(const double &aMix);
  void setIETMapping(const std::string &IETMapping);
  void setNIter(const int &nIter);
  void setOutIter(const int &outIter);
  void setRecoveryFileName(const std::string &recoveryFileName);
  void setGuess(const SlfcGuess &guess);
  // Getters
  double getErrMin() const { return errMin; }
  std::string getIETMapping() const { return IETMapping; }
  double getMixingParameter() const { return aMix; }
  int getNIter() const { return nIter; }
  int getOutIter() const { return outIter; }
  std::string getRecoveryFileName() const { return recoveryFileName; }
  SlfcGuess getGuess() const { return guess; }
  // Print content of the data structure
  void print() const;
  // Compare two StlsInput objects
  bool isEqual(const StlsInput &in) const;

protected:

  // Mixing parameter for the iterative procedure
  double aMix;
  // Minimum error for convergence in the iterative procedure
  double errMin;
  // Maximum number of iterations
  int nIter;
  // Output frequency
  int outIter;
  // Mapping between the quantum and classical state points for the IET-based
  // schemes
  std::string IETMapping;
  // Name of the file used to store the recovery data
  std::string recoveryFileName;
  // Initial guess
  SlfcGuess guess;
};

// -----------------------------------------------------------------
// Class to handle input for the QSTLS and QSTLS-IET schemes
// -----------------------------------------------------------------

class QstlsInput : public StlsInput {

public:

  // Typdef
  struct QstlsGuess {
    std::vector<double> wvg;
    std::vector<double> ssf;
    vecUtil::Vector2D adr;
    int matsubara = 0;
    bool operator==(const QstlsGuess &other) const {
      return wvg == other.wvg && ssf == other.ssf && adr == other.adr &&
             matsubara == other.matsubara;
    }
  };
  // Contructors
  explicit QstlsInput()
      : fixed(""),
        fixedIet("") {}
  // Setters
  void setFixed(const std::string &fixed);
  void setFixedIet(const std::string &fixedIet);
  void setGuess(const QstlsGuess &guess);
  // Getters
  std::string getFixed() const { return fixed; }
  std::string getFixedIet() const { return fixedIet; }
  QstlsGuess getGuess() const { return guess; }
  // Print content of the data structure
  void print() const;
  // Compare two QstlsInput objects
  bool isEqual(const QstlsInput &in) const;

private:

  // Name of the file with the fixed component of the auxiliary density response
  // (adr)
  std::string fixed;
  // Name of the file with the fixed component of the adr for iet schemes
  std::string fixedIet;
  // Initial guess
  QstlsGuess guess;
};

// -----------------------------------------------------------------
// Class to handle input for the VS schemes
// -----------------------------------------------------------------

class VSInput {

public:

  // Typdef
  struct FreeEnergyIntegrand {
    std::vector<double> grid;
    std::vector<double> alpha;
    std::vector<std::vector<double>> integrand;
    bool operator==(const FreeEnergyIntegrand &other) const {
      return grid == other.grid && integrand == other.integrand &&
             alpha == other.alpha;
    }
  };
  // Contructor
  explicit VSInput()
      : alphaGuess(std::vector<double>(2, 0)),
        drs(0),
        dTheta(0),
        errMinAlpha(0),
        nIterAlpha(0) {}
  // Setters
  void setAlphaGuess(const std::vector<double> &alphaGuess);
  void setCouplingResolution(const double &drs);
  void setDegeneracyResolution(const double &dTheta);
  void setErrMinAlpha(const double &errMinAlpha);
  void setNIterAlpha(const int &nIterAlpha);
  void setFreeEnergyIntegrand(const FreeEnergyIntegrand &freeEnergyIntegrand);
  // Getters
  std::vector<double> getAlphaGuess() const { return alphaGuess; }
  double getCouplingResolution() const { return drs; }
  double getDegeneracyResolution() const { return dTheta; }
  double getErrMinAlpha() const { return errMinAlpha; }
  double getNIterAlpha() const { return nIterAlpha; }
  FreeEnergyIntegrand getFreeEnergyIntegrand() const { return fxcIntegrand; }
  // Print content of the data structure
  void print() const;
  // Compare two VSStls objects
  bool isEqual(const VSInput &in) const;

private:

  // Initial guess for the free parameter
  std::vector<double> alphaGuess;
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
};

// -----------------------------------------------------------------
// Class to handle input for the VSStls scheme
// -----------------------------------------------------------------

class VSStlsInput : public StlsInput, public VSInput {

public:

  // Print content of the data structure
  void print() const;
  // Compare two VSStls objects
  bool isEqual(const VSStlsInput &in) const;
};

// -----------------------------------------------------------------
// Class to handle input for the QVSStls scheme
// -----------------------------------------------------------------

class QVSStlsInput : public QstlsInput, public VSInput {

public:

  // Print content of the data structure
  void print() const;
  // Compare two VSStls objects
  bool isEqual(const QVSStlsInput &in) const;
};

#endif
