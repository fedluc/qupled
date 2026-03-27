#ifndef INPUT_HPP
#define INPUT_HPP

#include "util/database.hpp"
#include "util/dimensions_util.hpp"
#include "util/num_util.hpp"
#include "util/vector2D.hpp"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------
// Default values
// -----------------------------------------------------------------

/** @brief Default value for uninitialized double parameters (signaling NaN). */
constexpr double DEFAULT_DOUBLE = numUtil::NaN;
/** @brief Default value for uninitialized integer parameters (signaling NaN).
 */
constexpr int DEFAULT_INT = numUtil::iNaN;
/** @brief Default value for uninitialized boolean parameters. */
constexpr bool DEFAULT_BOOL = false;

/**
 * @brief Base input class for all dielectric scheme solvers.
 *
 * Stores the thermodynamic state point (coupling parameter @p rs and
 * degeneracy parameter @p Theta), numerical grid parameters, and metadata
 * such as the theory name and database information. All derived input classes
 * extend this base with scheme-specific settings.
 */
class Input {

public:

  /** @brief Construct with all parameters set to their sentinel defaults. */
  explicit Input()
      : intError(DEFAULT_DOUBLE),
        rs(DEFAULT_DOUBLE),
        Theta(DEFAULT_DOUBLE),
        nThreads(DEFAULT_INT),
        isClassicTheory(DEFAULT_BOOL),
        isQuantumTheory(DEFAULT_BOOL),
        dimension(dimensionsUtil::Dimension::Default),
        dx(DEFAULT_DOUBLE),
        xmax(DEFAULT_DOUBLE),
        OmegaMax(DEFAULT_DOUBLE),
        nl(DEFAULT_INT) {}

  /** @brief Virtual destructor. */
  virtual ~Input() = default;

  /**
   * @brief Set the coupling parameter (Wigner–Seitz radius).
   * @param rs Non-negative coupling parameter.
   */
  void setCoupling(const double &rs);

  /**
   * @brief Set the database metadata for result storage.
   * @param dbInfo Database information struct.
   */
  void setDatabaseInfo(const databaseUtil::DatabaseInfo &dbInfo);

  /**
   * @brief Set the spatial dimension of the system.
   * @param dimension Dimension enum value (D2 or D3).
   */
  void setDimension(const dimensionsUtil::Dimension &dimension);

  /**
   * @brief Set the degeneracy parameter (reduced temperature @f$\Theta = k_B T
   * / E_F@f$).
   * @param Theta Non-negative degeneracy parameter.
   */
  void setDegeneracy(const double &Theta);

  /**
   * @brief Set the 2D integration quadrature scheme.
   * @param int2DScheme Name of the scheme (e.g., "segregated").
   */
  void setInt2DScheme(const std::string &int2DScheme);

  /**
   * @brief Set the relative error tolerance for numerical integrals.
   * @param intError Positive relative error tolerance.
   */
  void setIntError(const double &intError);

  /**
   * @brief Set the number of OpenMP threads for parallel sections.
   * @param nThreads Positive thread count.
   */
  void setNThreads(const int &nThreads);

  /**
   * @brief Set the name of the dielectric theory to solve.
   * @param theory Theory identifier string (e.g., "STLS", "qSTLS").
   */
  void setTheory(const std::string &theory);

  /**
   * @brief Set the initial bracket for the chemical potential root-finder.
   * @param muGuess Two-element vector [lower, upper] bracketing the solution.
   */
  void setChemicalPotentialGuess(const std::vector<double> &muGuess);

  /**
   * @brief Set the number of Matsubara frequencies used in finite-temperature
   * sums.
   * @param nMatsubara Positive integer number of frequencies.
   */
  void setNMatsubara(const int &nMatsubara);

  /**
   * @brief Set the wave-vector grid resolution.
   * @param dx Positive grid spacing.
   */
  void setWaveVectorGridRes(const double &dx);

  /**
   * @brief Set the wave-vector grid cutoff.
   * @param xmax Positive maximum wave-vector value.
   */
  void setWaveVectorGridCutoff(const double &xmax);

  /**
   * @brief Set the frequency cutoff (ground-state calculations only).
   * @param OmegaMax Positive maximum real frequency.
   */
  void setFrequencyCutoff(const double &OmegaMax);

  /** @brief Return the coupling parameter. */
  double getCoupling() const { return rs; }
  /** @brief Return the database metadata. */
  databaseUtil::DatabaseInfo getDatabaseInfo() const { return dbInfo; }
  /** @brief Return the spatial dimension. */
  dimensionsUtil::Dimension getDimension() const { return dimension; }
  /** @brief Return the degeneracy parameter. */
  double getDegeneracy() const { return Theta; }
  /** @brief Return the 2D integration quadrature scheme name. */
  std::string getInt2DScheme() const { return int2DScheme; }
  /** @brief Return the relative integral error tolerance. */
  double getIntError() const { return intError; }
  /** @brief Return the number of OpenMP threads. */
  int getNThreads() const { return nThreads; }
  /** @brief Return the theory name. */
  std::string getTheory() const { return theory; }
  /** @brief Return true if the theory is a classical (non-quantum)
   * approximation. */
  bool isClassic() const { return isClassicTheory; }
  /** @brief Return the chemical potential initial guess bracket. */
  std::vector<double> getChemicalPotentialGuess() const { return muGuess; }
  /** @brief Return the number of Matsubara frequencies. */
  int getNMatsubara() const { return nl; }
  /** @brief Return the wave-vector grid resolution. */
  double getWaveVectorGridRes() const { return dx; }
  /** @brief Return the wave-vector grid cutoff. */
  double getWaveVectorGridCutoff() const { return xmax; }
  /** @brief Return the frequency cutoff. */
  double getFrequencyCutoff() const { return OmegaMax; }

protected:

  /** @brief Relative error tolerance for numerical integrals. */
  double intError;
  /** @brief Quantum coupling parameter (Wigner–Seitz radius). */
  double rs;
  /** @brief Degeneracy parameter @f$\Theta = k_B T / E_F@f$. */
  double Theta;
  /** @brief Number of OpenMP threads for parallel sections. */
  int nThreads;
  /** @brief True if the theory is a classical (non-quantum) approximation. */
  bool isClassicTheory;
  /** @brief True if the theory includes quantum (exchange-correlation) effects.
   */
  bool isQuantumTheory;
  /** @brief 2D quadrature scheme name. */
  std::string int2DScheme;
  /** @brief Theory identifier string. */
  std::string theory;
  /** @brief Database metadata for result storage. */
  databaseUtil::DatabaseInfo dbInfo;
  /** @brief Spatial dimension (default: 3D). */
  dimensionsUtil::Dimension dimension;
  /** @brief Wave-vector grid resolution. */
  double dx;
  /** @brief Wave-vector grid cutoff. */
  double xmax;
  /** @brief Frequency cutoff (relevant only for ground-state calculations). */
  double OmegaMax;
  /** @brief Number of Matsubara frequencies. */
  int nl;
  /** @brief Initial bracket for the chemical potential root-finder. */
  std::vector<double> muGuess;
};

/**
 * @brief Carries the initial guess for iterative schemes.
 *
 * Holds a pre-computed wave-vector grid, static structure factor, and
 * local field correction that can be used to warm-start the iterative solver.
 */
struct Guess {
  /** @brief Wave-vector grid. */
  std::vector<double> wvg;
  /** @brief Static structure factor values. */
  std::vector<double> ssf;
  /** @brief Local field correction values (wave-vectors × Matsubara
   * frequencies). */
  Vector2D lfc;
};

/**
 * @brief Input class for iteratively solved dielectric schemes.
 *
 * Extends @p Input with convergence criteria (minimum error, maximum
 * iterations) and a linear-mixing parameter for the iterative updates.
 * An optional initial guess can be provided to warm-start the iteration.
 */
class IterationInput : public Input {

public:

  /** @brief Construct with all iteration parameters set to their sentinel
   * defaults. */
  explicit IterationInput()
      : aMix(DEFAULT_DOUBLE),
        errMin(DEFAULT_DOUBLE),
        nIter(DEFAULT_INT) {}

  /**
   * @brief Set the convergence threshold.
   * @param errMin Positive minimum RMS error for stopping the iteration.
   */
  void setErrMin(const double &errMin);

  /**
   * @brief Set the initial guess for the iterative solver.
   * @param guess Struct containing wvg, ssf, and lfc arrays.
   */
  void setGuess(const Guess &guess);

  /**
   * @brief Set the linear-mixing parameter.
   * @param aMix Mixing fraction in (0, 1].
   */
  void setMixingParameter(const double &aMix);

  /**
   * @brief Set the maximum number of iterations.
   * @param nIter Positive integer maximum.
   */
  void setNIter(const int &nIter);

  /** @brief Return the convergence threshold. */
  double getErrMin() const { return errMin; }
  /** @brief Return the initial guess struct. */
  Guess getGuess() const { return guess; }
  /** @brief Return the linear-mixing parameter. */
  double getMixingParameter() const { return aMix; }
  /** @brief Return the maximum number of iterations. */
  int getNIter() const { return nIter; }

protected:

  /** @brief Linear-mixing parameter for SSF updates. */
  double aMix;
  /** @brief Minimum RMS error for convergence. */
  double errMin;
  /** @brief Maximum number of self-consistency iterations. */
  int nIter;
  /** @brief Optional warm-start initial guess. */
  Guess guess;
};

/**
 * @brief Mixin input class for quantum (qSTLS and qSTLS-IET) schemes.
 *
 * Carries the run identifier used to load a pre-computed fixed component
 * of the auxiliary density response from the database.
 */
class QuantumInput {

public:

  /**
   * @brief Set the run ID for the fixed ADR component.
   * @param fixedRunId Integer run identifier in the database.
   */
  void setFixedRunId(const int &fixedRunId);

  /** @brief Return the fixed ADR run identifier. */
  int getFixedRunId() const { return fixedRunId; }

protected:

  /** @brief Database run ID for the fixed ADR component. */
  int fixedRunId;
};

/**
 * @brief Mixin input class for IET-based schemes.
 *
 * Carries the name of the quantum-to-classical mapping used when evaluating
 * the bridge function (e.g., "HNC", "IOI", "LCT").
 */
class IetInput {

public:

  /** @brief Virtual destructor. */
  virtual ~IetInput() = default;

  /**
   * @brief Set the quantum-classical mapping for the bridge function.
   * @param mapping Mapping identifier string.
   */
  void setMapping(const std::string &mapping);

  /** @brief Return the mapping identifier. */
  std::string getMapping() const { return mapping; }

protected:

  /** @brief Quantum-classical mapping identifier for the IET bridge function.
   */
  std::string mapping;
};

/**
 * @brief Input class for the STLS scheme.
 *
 * Inherits all iteration parameters from @p IterationInput with no
 * additional STLS-specific settings.
 */
class StlsInput : public IterationInput {

public:

  /** @brief Construct with default parameters. */
  explicit StlsInput() = default;
};

/**
 * @brief Input class for the STLS-IET scheme.
 *
 * Combines STLS iteration parameters with the IET mapping setting.
 */
class StlsIetInput : public StlsInput, public IetInput {

public:

  /** @brief Construct with default parameters. */
  explicit StlsIetInput() = default;
};

/**
 * @brief Input class for the qSTLS scheme.
 *
 * Combines STLS iteration parameters with the fixed ADR run identifier.
 */
class QstlsInput : public StlsInput, public QuantumInput {

public:

  /** @brief Construct with default parameters. */
  explicit QstlsInput() = default;
};

/**
 * @brief Input class for the qSTLS-IET scheme.
 *
 * Combines qSTLS parameters with the IET mapping setting.
 */
class QstlsIetInput : public QstlsInput, public IetInput {

public:

  /** @brief Construct with default parameters. */
  explicit QstlsIetInput() = default;
};

/**
 * @brief Mixin input class for variational-swarm (VS) schemes.
 *
 * Carries the free-parameter guess, finite-difference grid resolutions,
 * convergence settings for the alpha iterations, and optionally a
 * pre-computed free-energy integrand to warm-start the VS optimization.
 */
class VSInput {

public:

  /**
   * @brief Pre-computed free-energy integrand for warm-starting.
   *
   * The @p grid field is the coupling-parameter grid and @p integrand
   * contains the per-theta integrand values.
   */
  struct FreeEnergyIntegrand {
    /** @brief Coupling-parameter grid. */
    std::vector<double> grid;
    /** @brief Free-energy integrand values (indexed by theta, then rs). */
    std::vector<std::vector<double>> integrand;
  };

  /** @brief Construct with all VS parameters set to their sentinel defaults. */
  explicit VSInput()
      : drs(DEFAULT_DOUBLE),
        dTheta(DEFAULT_DOUBLE),
        errMinAlpha(DEFAULT_DOUBLE),
        nIterAlpha(DEFAULT_INT) {}

  /** @brief Virtual destructor. */
  virtual ~VSInput() = default;

  /**
   * @brief Set the initial bracket for the free parameter alpha.
   * @param alphaGuess Two-element vector bracketing alpha.
   */
  void setAlphaGuess(const std::vector<double> &alphaGuess);

  /**
   * @brief Set the coupling-parameter grid resolution for finite differences.
   * @param drs Positive grid spacing in rs.
   */
  void setCouplingResolution(const double &drs);

  /**
   * @brief Set the degeneracy-parameter grid resolution for finite differences.
   * @param dTheta Positive grid spacing in Theta.
   */
  void setDegeneracyResolution(const double &dTheta);

  /**
   * @brief Set the convergence threshold for the alpha iterations.
   * @param errMinAlpha Positive minimum error for stopping.
   */
  void setErrMinAlpha(const double &errMinAlpha);

  /**
   * @brief Set the maximum number of alpha iterations.
   * @param nIterAlpha Positive integer maximum.
   */
  void setNIterAlpha(const int &nIterAlpha);

  /**
   * @brief Provide a pre-computed free-energy integrand for warm-starting.
   * @param freeEnergyIntegrand Struct containing grid and integrand values.
   */
  void setFreeEnergyIntegrand(const FreeEnergyIntegrand &freeEnergyIntegrand);

  /** @brief Return the alpha guess bracket. */
  std::vector<double> getAlphaGuess() const { return alphaGuess; }
  /** @brief Return the coupling-parameter grid resolution. */
  double getCouplingResolution() const { return drs; }
  /** @brief Return the degeneracy-parameter grid resolution. */
  double getDegeneracyResolution() const { return dTheta; }
  /** @brief Return the alpha convergence threshold. */
  double getErrMinAlpha() const { return errMinAlpha; }
  /** @brief Return the maximum number of alpha iterations. */
  double getNIterAlpha() const { return nIterAlpha; }
  /** @brief Return the pre-computed free-energy integrand. */
  FreeEnergyIntegrand getFreeEnergyIntegrand() const { return fxcIntegrand; }

private:

  /** @brief Initial bracket for the free parameter alpha. */
  std::vector<double> alphaGuess;
  /** @brief Coupling-parameter grid resolution (finite-difference step in rs).
   */
  double drs;
  /** @brief Degeneracy-parameter grid resolution (finite-difference step in
   * Theta). */
  double dTheta;
  /** @brief Convergence threshold for the alpha iterations. */
  double errMinAlpha;
  /** @brief Maximum number of alpha iterations. */
  int nIterAlpha;
  /** @brief Pre-computed free-energy integrand for warm-starting. */
  FreeEnergyIntegrand fxcIntegrand;
};

/**
 * @brief Input class for the VS-STLS scheme.
 *
 * Combines VS optimization parameters with STLS iteration parameters.
 */
class VSStlsInput : public VSInput, public StlsInput {

public:

  /** @brief Construct with default parameters. */
  explicit VSStlsInput() = default;
};

/**
 * @brief Input class for the QVS-STLS (quantum VS-STLS) scheme.
 *
 * Combines VS optimization parameters with qSTLS iteration parameters
 * (including the fixed ADR run identifier).
 */
class QVSStlsInput : public VSInput, public QstlsInput {

public:

  /** @brief Construct with default parameters. */
  explicit QVSStlsInput() = default;
};

#endif
