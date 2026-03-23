#ifndef VS_VSBASE_HPP
#define VS_VSBASE_HPP

#include "input.hpp"
#include "logger.hpp"
#include "numerics.hpp"
#include "vs/grid_point.hpp"
#include <vector>

// Forward declaration
class VSManager;

// -----------------------------------------------------------------
// QAdder: computes the Q-term in the VS free parameter expression
// -----------------------------------------------------------------

class QAdder {
public:

  enum class Mode { CLASSICAL, QUANTUM };

  static QAdder classical(const std::vector<double> &wvg,
                          const std::vector<double> &ssf,
                          const std::shared_ptr<const Input> &in);

  static QAdder quantum(double Theta,
                        double mu,
                        double limitMin,
                        double limitMax,
                        const std::vector<double> &itgGrid,
                        std::shared_ptr<Integrator1D> itg1,
                        std::shared_ptr<Integrator2D> itg2,
                        std::shared_ptr<Interpolator1D> interp);

  double get() const;

private:

  Mode mode;

  // Classical fields
  std::vector<double> classWvg;
  std::vector<double> classSsf;
  std::shared_ptr<const Input> classIn;

  // Quantum fields
  double Theta;
  double mu;
  std::pair<double, double> limits;
  const std::vector<double> *itgGridPtr;
  std::shared_ptr<Integrator1D> itg1;
  std::shared_ptr<Integrator2D> itg2;
  std::shared_ptr<Interpolator1D> interp;

  // Quantum-only helpers
  double ssf(const double &y) const;
  double integrandDenominator(const double y) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  void getIntDenominator(double &res) const;
};

// -----------------------------------------------------------------
// VSBase class
// -----------------------------------------------------------------

class VSBase : public Logger {

public:

  explicit VSBase() {}
  virtual ~VSBase() = default;
  // Solve the VS scheme
  int compute();
  // Getters
  double getAlpha() const { return alpha; }
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const;
  const std::vector<double> &getFreeEnergyGrid() const;
  // Public output interface (delegates to grid())
  const std::vector<double> &getSsf() const;
  const Vector2D &getLfc() const;
  const std::vector<double> &getWvg() const;
  const Vector2D &getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  double getError() const;

protected:

  // Free parameter
  double alpha;
  // Free energy integrand data (indexed by theta, then rs)
  std::vector<double> rsGrid;
  // Free energy integrand values at the 3 theta grid points (DOWN, CENTER, UP)
  std::vector<std::vector<double>> fxcIntegrand;
  // Abstract interface
  virtual const VSInput &in() const = 0;
  virtual const Input &inScheme() const = 0;
  // Virtual grid accessor
  virtual VSManager &grid() = 0;
  virtual const VSManager &grid() const = 0;
  // Grid coordination (non-virtual, uses grid())
  int runGrid();
  double getCoupling(GridPoint p) const;
  double getDegeneracy(GridPoint p) const;
  double getFxcIntegrandValue(GridPoint p) const;
  // Scheme-specific Q computation (virtual, overridden by derived classes)
  virtual double computeQRaw(GridPoint p) const = 0;
  // Non-virtual helpers
  void setRsGrid();
  void setFxcIntegrand();
  void updateFxcIntegrand();
  double computeFreeEnergy(GridPoint p, bool normalize) const;
  std::vector<double> getFreeEnergyData() const;
  GridPoint getOutputGridPoint() const;

private:

  double computeAlpha();
  std::vector<double> computeQData();
  void doIterations();
  double alphaDifference(const double &alphaTmp);
};

#endif
