#ifndef VSSTLS_NEW_HPP
#define VSSTLS_NEW_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vsbase_new.hpp"
#include <limits>
#include <map>

class ThermoProp: public ThermoPropBase {

public:
  
  // Constructor
  explicit ThermoProp(const VSStlsInput &in): ThermoPropBase(in) {}
  
};


class VSStlsNew : public VSBase, public Rpa {

public:

  // Constructor from initial data
  explicit VSStlsNew(const VSStlsInput &in_)
    : VSBase(in_), Rpa(in_), thermoProp(in_) {}
  // Constructor for recursive calculations
  VSStlsNew(const VSStlsInput &in_, const ThermoProp &thermoProp_)
    : VSBase(in_), Rpa(in_, false), thermoProp(in_) {
    thermoProp.copyFreeEnergyIntegrand(thermoProp_);
  }
  
  // Solve the scheme
  using VSBase::compute;
  
private:

  // Input
  VSStlsInput in;
  // Verbosity
  using VSBase::verbose;
  // Getter
  const ThermoPropBase& getThermoProp() const { return thermoProp; }
  ThermoPropBase& getThermoProp() { return thermoProp; }
  // Thermodynamic properties
  ThermoProp thermoProp;
  // Initialize
  void initScheme();
  void initFreeEnergyIntegrand();
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
};

#endif
