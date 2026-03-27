#ifndef VS_VSBASE_HPP
#define VS_VSBASE_HPP

#include "schemes/input.hpp"
#include "util/logger.hpp"
#include "util/numerics.hpp"
#include "vs/grid_point.hpp"
#include <vector>

// Forward declaration
class VSManager;

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
  virtual VSManager &grid() = 0;
  virtual const VSManager &grid() const = 0;
  // Grid coordination (non-virtual, uses grid())
  int runGrid();
  double getCoupling(GridPoint p) const;
  double getDegeneracy(GridPoint p) const;
  double getFxcIntegrandValue(GridPoint p) const;
  double getQAdder(GridPoint p) const;
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
