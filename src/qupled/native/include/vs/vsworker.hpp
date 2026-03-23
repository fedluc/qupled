#ifndef VS_VSWORKER_HPP
#define VS_VSWORKER_HPP

#include "vector2D.hpp"
#include <vector>

// -----------------------------------------------------------------
// VSWorker: abstract interface for all worker objects
// -----------------------------------------------------------------

class VSWorker {
public:

  virtual ~VSWorker() = default;
  // Getters
  virtual const Vector2D &getLfc() const = 0;
  virtual const std::vector<double> &getWvg() const = 0;
  virtual const std::vector<double> &getSsf() const = 0;
  // Iteration protocol
  virtual void init() = 0;
  virtual void initialGuess() = 0;
  virtual void computeSsf() = 0;
  // Synchronized LFC protocol
  virtual void computeLfc() = 0;
  virtual void applyLfcDiff(const Vector2D &v) = 0;
  virtual double computeError() const = 0;
  virtual void updateSolution() = 0;
};

#endif
