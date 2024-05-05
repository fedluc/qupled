#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

class ChemicalPotential {

public:

  explicit ChemicalPotential(const double &Theta_)
      : Theta(Theta_) {}
  void compute(const std::vector<double> &guess);
  double get() const { return mu; }

private:

  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  double mu;
  // Normalization condition
  double normalizationCondition(const double &mu) const;
};

#endif
