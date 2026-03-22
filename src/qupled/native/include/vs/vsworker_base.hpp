#ifndef VS_VSWORKER_BASE_HPP
#define VS_VSWORKER_BASE_HPP

#include "vector2D.hpp"
#include <vector>

// -----------------------------------------------------------------
// VSWorkerBase: abstract interface for all worker objects
// -----------------------------------------------------------------

class VSWorkerBase {
public:

  virtual ~VSWorkerBase() = default;
  // Synchronized LFC protocol
  virtual void computeBaseLfc() = 0;
  virtual void applyLfcDiff(const Vector2D &v) = 0;
  // Getters
  virtual const Vector2D &getLfc() const = 0;
  virtual const std::vector<double> &getWvg() const = 0;
  virtual const std::vector<double> &getSsf() const = 0;
  // Iteration protocol
  virtual void init() = 0;
  virtual void initialGuess() = 0;
  virtual void computeSsf() = 0;
  virtual double computeError() const = 0;
  virtual void updateSolution() = 0;
};

#endif
