#ifndef VSBASE_HPP
#define VSBASE_HPP

#include "input.hpp"
#include "logger.hpp"
#include "numerics.hpp"
#include "state_point_grid.hpp"
#include <vector>

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
  virtual int runGrid() = 0;
  virtual double getCoupling(GridPoint p) const = 0;
  virtual double getDegeneracy(GridPoint p) const = 0;
  virtual double getFxcIntegrandValue(GridPoint p) const = 0;
  virtual std::vector<double> computeQData() = 0;
  // Non-virtual helpers
  void setRsGrid();
  void setFxcIntegrand();
  void updateFxcIntegrand();
  double computeFreeEnergy(GridPoint p, bool normalize) const;
  std::vector<double> getFreeEnergyData() const;
  GridPoint getOutputGridPoint() const;

private:

  double computeAlpha();
  void doIterations();
  double alphaDifference(const double &alphaTmp);
};

#endif
