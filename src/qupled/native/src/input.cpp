#include "input.hpp"
#include "mpi_util.hpp"
#include <cmath>
#include <fmt/core.h>

using namespace std;
using namespace MPIUtil;

// -----------------------------------------------------------------
// Input class
// -----------------------------------------------------------------

void Input::setCoupling(const double &rs_) {
  if (rs_ < 0) {
    throwError("The quantum coupling parameter can't be negative");
  }
  this->rs = rs_;
}

void Input::setDatabaseName(const string &databaseName) {
  this->databaseName = databaseName;
}

void Input::setDatabaseRunId(const int &databaseRunId) {
  this->databaseRunId = databaseRunId;
}

void Input::setDegeneracy(const double &Theta_) {
  if (Theta_ < 0.0) {
    throwError("The quantum degeneracy parameter can't be negative");
  }
  this->Theta = Theta_;
}

void Input::setTheory(const string &theory_) {
  const vector<string> cTheories = {
      "RPA", "ESA", "STLS", "STLS-HNC", "STLS-IOI", "STLS-LCT", "VSSTLS"};
  const vector<string> qTheories = {
      "QSTLS", "QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT", "QVSSTLS"};
  isClassicTheory = count(cTheories.begin(), cTheories.end(), theory_) != 0;
  isQuantumTheory = count(qTheories.begin(), qTheories.end(), theory_) != 0;
  if (!isClassicTheory && !isQuantumTheory) {
    throwError("Invalid dielectric theory: " + theory_);
  }
  // A theory can't both be classical and quantum at the same time
  assert(!isClassicTheory || !isQuantumTheory);
  this->theory = theory_;
}

void Input::setInt2DScheme(const string &int2DScheme) {
  const vector<string> schemes = {"full", "segregated"};
  if (count(schemes.begin(), schemes.end(), int2DScheme) == 0) {
    throwError("Unknown scheme for 2D integrals: " + int2DScheme);
  }
  this->int2DScheme = int2DScheme;
}

void Input::setIntError(const double &intError) {
  if (intError <= 0) {
    throwError(
        "The accuracy for the integral computations must be larger than zero");
  }
  this->intError = intError;
}

void Input::setNThreads(const int &nThreads) {
  if (nThreads <= 0) {
    throwError("The number of threads must be larger than zero");
  }
  this->nThreads = nThreads;
}

// -----------------------------------------------------------------
// IterationInput class
// -----------------------------------------------------------------

void IterationInput::setErrMin(const double &errMin) {
  if (errMin <= 0.0) {
    throwError("The minimum error for convergence must be larger than zero");
  }
  this->errMin = errMin;
}

void IterationInput::setMixingParameter(const double &aMix) {
  if (aMix < 0.0 || aMix > 1.0) {
    throwError("The mixing parameter must be a number between zero and one");
  }
  this->aMix = aMix;
}

void IterationInput::setNIter(const int &nIter) {
  if (nIter < 0) {
    throwError("The maximum number of iterations can't be negative");
  }
  this->nIter = nIter;
}

void IterationInput::setOutIter(const int &outIter) {
  if (outIter < 0) { throwError("The output frequency can't be negative"); }
  this->outIter = outIter;
}

// -----------------------------------------------------------------
// ClassicInput class
// -----------------------------------------------------------------

void ClassicInput::setIETMapping(const string &IETMapping) {
  const vector<string> mappings = {"standard", "sqrt", "linear"};
  if (count(mappings.begin(), mappings.end(), IETMapping) == 0) {
    throwError("Unknown IET mapping: " + IETMapping);
  }
  this->IETMapping = IETMapping;
}

void ClassicInput::setGuess(const Guess &guess) {
  if (guess.wvg.size() != guess.slfc.size()) {
    throwError("The initial guess is inconsistent");
  }
  this->guess = guess;
}

// -----------------------------------------------------------------
// QuantumInput class
// -----------------------------------------------------------------

void QuantumInput::setFixed(const string &fixed) { this->fixed = fixed; }

void QuantumInput::setFixedIet(const string &fixedIet) {
  this->fixedIet = fixedIet;
}

void QuantumInput::setGuess(const Guess &guess) {
  const size_t nl = guess.matsubara;
  const size_t nx = guess.wvg.size();
  const size_t adrRows = (guess.adr.empty()) ? 0 : nx;
  const size_t adrCols = (guess.adr.empty()) ? 0 : nl;
  const bool consistentGuess = guess.ssf.size() == nx
                               && guess.adr.size(0) == adrRows
                               && guess.adr.size(1) == adrCols;
  if (!consistentGuess) { throwError("The initial guess is inconsistent"); }
  this->guess = guess;
}

// -----------------------------------------------------------------
// RpaInput class
// -----------------------------------------------------------------

void RpaInput::setChemicalPotentialGuess(const vector<double> &muGuess) {
  if (muGuess.size() != 2 || muGuess[0] >= muGuess[1]) {
    throwError("Invalid guess for chemical potential calculation");
  }
  this->muGuess = muGuess;
}

void RpaInput::setNMatsubara(const int &nl) {
  if (nl < 0) {
    throwError("The number of matsubara frequencies can't be negative");
  }
  this->nl = nl;
}

void RpaInput::setWaveVectorGridRes(const double &dx) {
  if (dx <= 0.0) {
    throwError("The wave-vector grid resolution must be larger than zero");
  }
  this->dx = dx;
}

void RpaInput::setWaveVectorGridCutoff(const double &xmax) {
  if (xmax <= 0.0) {
    throwError("The wave-vector grid cutoff must be larger than zero");
  }
  this->xmax = xmax;
}

void RpaInput::setFrequencyCutoff(const double &OmegaMax) {
  if (OmegaMax <= 0.0) {
    throwError("The frequency cutoff must be larger than zero");
  }
  this->OmegaMax = OmegaMax;
}

// -----------------------------------------------------------------
// VSInput class
// -----------------------------------------------------------------

void VSInput::setCouplingResolution(const double &drs) {
  if (drs <= 0) {
    throwError("The coupling parameter resolution must be larger than zero");
  }
  this->drs = drs;
}

void VSInput::setDegeneracyResolution(const double &dTheta) {
  if (dTheta <= 0) {
    throwError("The degeneracy parameter resolution must be larger than zero");
  }
  this->dTheta = dTheta;
}

void VSInput::setAlphaGuess(const vector<double> &alphaGuess) {
  if (alphaGuess.size() != 2 || alphaGuess[0] >= alphaGuess[1]) {
    throwError("Invalid guess for free parameter calculation");
  }
  this->alphaGuess = alphaGuess;
}

void VSInput::setErrMinAlpha(const double &errMinAlpha) {
  if (errMinAlpha <= 0.0) {
    throwError("The minimum error for convergence must be larger than zero");
  }
  this->errMinAlpha = errMinAlpha;
}

void VSInput::setNIterAlpha(const int &nIterAlpha) {
  if (nIterAlpha < 0) {
    throwError("The maximum number of iterations can't be negative");
  }
  this->nIterAlpha = nIterAlpha;
}

void VSInput::setFreeEnergyIntegrand(const FreeEnergyIntegrand &fxcIntegrand) {
  const auto &integrands = fxcIntegrand.integrand;
  const size_t referenceSize = (integrands.empty()) ? 0 : integrands[0].size();
  for (const auto &integrand : integrands) {
    if (integrand.size() != referenceSize) {
      throwError("The free energy integrand is inconsistent");
    }
  }
  if (fxcIntegrand.grid.size() != referenceSize) {
    throwError("The free energy integrand is inconsistent");
  }
  this->fxcIntegrand = fxcIntegrand;
}