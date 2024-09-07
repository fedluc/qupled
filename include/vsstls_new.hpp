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
  
};


class ThermoProp : public ThermoPropBase {

public:

  // Constructor
  explicit ThermoProp(const VSStlsInput &in_);

private:

  // Structural properties
  std::shared_ptr<StructPropNew> structProp;
  
};


class VSStlsNew : public VSBase, public Stls {

public:

  // Constructor from initial data
  explicit VSStlsNew(const VSStlsInput &in_);
  // Constructor for recursive calculations
  VSStlsNew(const VSStlsInput &in_, const ThermoProp &thermoProp_);
  
  // Solve the scheme
  using VSBase::compute;
  
private:

  // Input
  VSStlsInput in;
  // Verbosity
  using VSBase::verbose;
  // Thermodynamic properties
  std::shared_ptr<ThermoProp> thermoProp;
  // Initialize
  void initScheme();
  void initFreeEnergyIntegrand();
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
};

  
#endif
