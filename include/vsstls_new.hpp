#ifndef VSSTLS_NEW_HPP
#define VSSTLS_NEW_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vsbase_new.hpp"
#include <limits>
#include <map>

class ThermoProp : public ThermoPropBase {
};

class VSStlsNew : public VSBase, public Rpa {

public:

  // Constructor from initial data
  explicit VSStlsNew(const VSStlsInput &in_)
    : VSStlsNew(in_, std::make_shared<ThermoProp>()) {}
  // Constructor for recursive calculations
  VSStlsNew(const VSStlsInput &in_, const std::shared_ptr<ThermoProp> &thermoProp_)
    : VSBase(in_, thermoProp_, false), Rpa(in_, false) {}

private:
  
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
};

#endif
