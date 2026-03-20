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

  double alpha;
  std::vector<double> rsGrid;
  // fxcIntegrand[0/1/2][rs] = free energy integrand at theta-down/center/up
  std::vector<std::vector<double>> fxcIntegrand;

  // Abstract interface implemented by VSStls / QVSStls
  virtual const VSInput &in() const = 0; // VS-specific input (alpha, drs, ...)
  virtual const Input &
  inScheme() const = 0; // base scheme input (rs, theta, ...)
  virtual void init() = 0;
  virtual void updateSolution() = 0;

  // Grid access — implemented by VSStls / QVSStls via their StatePointGrid
  virtual int runGrid() = 0;
  virtual double getCoupling(GridPoint p) const = 0;
  virtual double getDegeneracy(GridPoint p) const = 0;
  virtual double getFxcIntegrandValue(GridPoint p) const = 0;
  // Returns {Q, Qr, Qt}:
  //   VSStls  → {uint, uintr, uintt}  (internal energy and derivatives)
  //   QVSStls → quantum Q-term and derivatives
  virtual std::vector<double> computeQData() = 0;

  // Non-virtual helpers (logic previously in ThermoPropBase)
  void setRsGrid();
  void setFxcIntegrand();
  void updateFxcIntegrand();
  double computeFreeEnergy(GridPoint p, bool normalize) const;
  std::vector<double> getFreeEnergyData() const;
  GridPoint getOutputGridPoint() const;

private:

  // Unified alpha formula — not virtual
  double computeAlpha();
  void doIterations();
  double alphaDifference(const double &alphaTmp);
};

#endif
