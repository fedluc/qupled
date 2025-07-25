#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include "dimensions_util.hpp"
#include "input.hpp"
#include <memory>

class ChemicalPotential : public dimensionsUtil::DimensionsHandler {
public:

  explicit ChemicalPotential(const std::shared_ptr<const Input> in_)
      : in(in_) {}
  double get() const { return mu; }

private:

  const std::shared_ptr<const Input> in;
  double mu = DEFAULT_DOUBLE;
  void compute2D() override;
  void compute3D() override;
  double normalizationCondition(const double &mu) const;
};

#endif