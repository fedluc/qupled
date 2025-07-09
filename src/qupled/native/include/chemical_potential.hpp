#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include "input.hpp"
#include "dimensions_util.hpp"
#include <memory>
#include <vector>

class ChemicalPotential : public dimensionsUtil::DimensionAware<ChemicalPotential> {
public:
    explicit ChemicalPotential(const std::shared_ptr<const Input> in_)
        : in(in_) {}
    
    double get() const { return mu; }
    
    void compute(Input::Dimension dim) {
        DimensionAware<ChemicalPotential>::compute(dim);
    }

private:
    const std::shared_ptr<const Input> in;
    double mu;
    
    friend class dimensionsUtil::DimensionAware<ChemicalPotential>;
    void compute2D();
    void compute3D();
    double normalizationCondition(const double &mu) const;
};

#endif