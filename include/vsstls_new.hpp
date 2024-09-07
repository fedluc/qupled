#ifndef VSSTLS_NEW_HPP
#define VSSTLS_NEW_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vsbase_new.hpp"
#include <limits>
#include <map>

class StlsCSRNew: public CSR, public Stls {

public:
  
  // Constructor
  explicit StlsCSRNew(const VSStlsInput &in_)
    : CSR(in_), Stls(in_.toStlsInput(), false, false) {}
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();

private:

  
};

class StructPropNew : public StructPropBase {

public:
  
  explicit StructPropNew(const VSStlsInput &in_) : StructPropBase(in_) {}

private:

  void doIterations();
  std::vector<CSR> csr;
  const std::vector<CSR>& getCsr() const { return csr; }
  std::vector<CSR>& getCsr() { return csr; }
  
  
};


class ThermoProp : public ThermoPropBase {

public:
  
  explicit ThermoProp(const VSStlsInput &in_): ThermoPropBase(in_), structProp(in_) {}

private:

  // Structural properties
  StructPropNew structProp;
  // Getter
  const StructPropBase& getStructProp() const { return structProp; }
  StructPropBase& getStructProp() { return structProp; }
  
};


class VSStlsNew : public VSBase, public Stls {

public:

  // Constructor from initial data
  explicit VSStlsNew(const VSStlsInput &in_)
    : VSBase(in_), Stls(in_.toStlsInput()), thermoProp(in_) {}
  // Constructor for recursive calculations
  VSStlsNew(const VSStlsInput &in_, const ThermoProp &thermoProp_)
    : VSBase(in_), Stls(in_.toStlsInput(), false, false), thermoProp(in_) {
    thermoProp.copyFreeEnergyIntegrand(thermoProp_);
  }
  
  // Solve the scheme
  using VSBase::compute;
  
private:

  // Input
  VSStlsInput in;
  // Verbosity
  using VSBase::verbose;
  // Thermodynamic properties
  ThermoProp thermoProp;
  // Getter
  const ThermoPropBase& getThermoProp() const { return thermoProp; }
  ThermoPropBase& getThermoProp() { return thermoProp; }
  // Initialize
  void initScheme();
  void initFreeEnergyIntegrand();
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
};

  
#endif
