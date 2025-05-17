#ifndef INPUT_HPP
#define INPUT_HPP

#include "database.hpp"
#include "num_util.hpp"
#include "vector2D.hpp"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------
// Default values
// -----------------------------------------------------------------

constexpr double DEFAULT_DOUBLE = numUtil::NaN;
constexpr int DEFAULT_INT = numUtil::iNaN;
constexpr bool DEFAULT_BOOL = false;

// -----------------------------------------------------------------
// Base class to handle input for the dielectric schemes
// -----------------------------------------------------------------
class Input {

public:

  // Constructor
  explicit Input()
      : intError(DEFAULT_DOUBLE),
        rs(DEFAULT_DOUBLE),
        Theta(DEFAULT_DOUBLE),
        nThreads(DEFAULT_INT),
        isClassicTheory(DEFAULT_BOOL),
        isQuantumTheory(DEFAULT_BOOL),
        dx(DEFAULT_DOUBLE),
        xmax(DEFAULT_DOUBLE),
        OmegaMax(DEFAULT_DOUBLE),
        nl(DEFAULT_INT) {}
  // Destructor
  virtual ~Input() = default;
  // Setters
  void setCoupling(const double &rs);
  void setDatabaseInfo(const DatabaseInfo &dbInfo);
  void setDegeneracy(const double &Theta);
  void setInt2DScheme(const std::string &int2DScheme);
  void setIntError(const double &intError);
  void setNThreads(const int &nThreads);
  void setTheory(const std::string &theory);
  void setChemicalPotentialGuess(const std::vector<double> &muGuess);
  void setNMatsubara(const int &nMatsubara);
  void setWaveVectorGridRes(const double &dx);
  void setWaveVectorGridCutoff(const double &xmax);
  void setFrequencyCutoff(const double &OmegaMax);

  // Getters
  double getCoupling() const { return rs; }
  DatabaseInfo getDatabaseInfo() const { return dbInfo; }
  double getDegeneracy() const { return Theta; }
  std::string getInt2DScheme() const { return int2DScheme; }
  double getIntError() const { return intError; }
  int getNThreads() const { return nThreads; }
  std::string getTheory() const { return theory; }
  bool isClassic() const { return isClassicTheory; }
  std::vector<double> getChemicalPotentialGuess() const { return muGuess; }
  int getNMatsubara() const { return nl; }
  double getWaveVectorGridRes() const { return dx; }
  double getWaveVectorGridCutoff() const { return xmax; }
  double getFrequencyCutoff() const { return OmegaMax; }

protected:

  // Accuracy for the integrals
  double intError;
  // Quantum coupling parameter
  double rs;
  // Degeneracy parameter
  double Theta;
  // Number of threads for parallel calculations
  int nThreads;
  // Type of theory
  bool isClassicTheory;
  bool isQuantumTheory;
  // Scheme for 2D integrals
  std::string int2DScheme;
  // Theory to be solved
  std::string theory;
  // Database information
  DatabaseInfo dbInfo;
  // Wave-vector grid resolution
  double dx;
  // Cutoff for the wave-vector grid
  double xmax;
  // Cutoff for the frequency (only relevant in the ground state)
  double OmegaMax;
  // Number of Matsubara frequencies
  int nl;
  // Initial guess for the chemical potential calculation
  std::vector<double> muGuess;
};

// -----------------------------------------------------------------
// Class to handle input for the schemes that are solved iteratively
// -----------------------------------------------------------------

// Typedef
struct Guess {
  std::vector<double> wvg;
  std::vector<double> ssf;
  Vector2D lfc;
};

class IterationInput : public Input {

public:

  // Contructor
  explicit IterationInput()
      : aMix(DEFAULT_DOUBLE),
        errMin(DEFAULT_DOUBLE),
        nIter(DEFAULT_INT) {}
  // Setters
  void setErrMin(const double &errMin);
  void setGuess(const Guess &guess);
  void setMixingParameter(const double &aMix);
  void setNIter(const int &nIter);
  // Getters
  double getErrMin() const { return errMin; }
  Guess getGuess() const { return guess; }
  double getMixingParameter() const { return aMix; }
  int getNIter() const { return nIter; }

protected:

  // Mixing parameter for the iterative procedure
  double aMix;
  // Minimum error for convergence in the iterative procedure
  double errMin;
  // Maximum number of iterations
  int nIter;
  // Initial guess
  Guess guess;
};

// -----------------------------------------------------------------
// Class to handle input for the QSTLS and QSTLS-IET schemes
// -----------------------------------------------------------------

class QuantumInput {

public:

  // Setters
  void setFixedRunId(const int &fixedRunId);
  void setFixedIet(const std::string &fixedIet);
  // Getters
  int getFixedRunId() const { return fixedRunId; }
  std::string getFixedIet() const { return fixedIet; }

protected:

  // Name of the file with the fixed component of the auxiliary density response
  // (adr)
  int fixedRunId;
  // Name of the file with the fixed component of the adr for iet schemes
  std::string fixedIet;
};

// -----------------------------------------------------------------
// Class to handle input for the iet schemes
// -----------------------------------------------------------------

class IetInput {

public:

  void setMapping(const std::string &mapping);
  // Getters
  std::string getMapping() const { return mapping; }

protected:

  // Mapping between the quantum and classical state points for the IET-based
  // schemes
  std::string mapping;
};

// -----------------------------------------------------------------
// Class to handle input for the STLS and STLS-IET schemes
// -----------------------------------------------------------------

class StlsInput : public IterationInput {

public:

  // Constructors
  explicit StlsInput() = default;
};

// -----------------------------------------------------------------
// Class to handle input for the STLS-IET schemes
// -----------------------------------------------------------------

class StlsIetInput : public StlsInput, public IetInput {

public:

  // Constructors
  explicit StlsIetInput() = default;
};

// -----------------------------------------------------------------
// Class to handle input for the QSTLS and QSTLS-IET schemes
// -----------------------------------------------------------------

class QstlsInput : public StlsInput, public QuantumInput {

public:

  // Constructors
  explicit QstlsInput() = default;
};

// -----------------------------------------------------------------
// Class to handle input for the QSTLS-IET schemes
// -----------------------------------------------------------------

class QstlsIetInput : public QstlsInput, public IetInput {

public:

  // Constructors
  explicit QstlsIetInput() = default;
};

// -----------------------------------------------------------------
// Class to handle input for the VS schemes
// -----------------------------------------------------------------

class VSInput {

public:

  // Typdef
  struct FreeEnergyIntegrand {
    std::vector<double> grid;
    std::vector<std::vector<double>> integrand;
  };
  // Contructor
  explicit VSInput()
      : drs(DEFAULT_DOUBLE),
        dTheta(DEFAULT_DOUBLE),
        errMinAlpha(DEFAULT_DOUBLE),
        nIterAlpha(DEFAULT_INT) {}
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

class VSStlsInput : public VSInput, public StlsInput {

public:

  // Constructors
  explicit VSStlsInput() = default;
};

// -----------------------------------------------------------------
// Class to handle input for the QVSStls scheme
// -----------------------------------------------------------------

class QVSStlsInput : public VSInput, public QstlsInput {

public:

  // Constructors
  explicit QVSStlsInput() = default;
};

#endif
