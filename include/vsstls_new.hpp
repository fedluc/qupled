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
    : CSR(in_), Stls(in_.toStlsInput(), false, false), in(in_) {}
  
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();

  // Publicly esposed private stls methods
  void init() { Stls::init(); }
  void initialGuess() { Stls::initialGuess(); }
  void computeSsf() { Stls::computeSsf(); }
  double computeError() { return Stls::computeError(); }
  void updateSolution() { Stls::updateSolution(); }

  // Getters
  std::vector<double> getSsf() const  { return Stls::getSsf(); }
  std::vector<double> getSlfc() const { return Stls::getSlfc(); }
  std::vector<double> getWvg() const { return Stls::getWvg(); }
  
private:

  // Input parameters
  VSStlsInput in;
  // Helper methods to compute the derivatives
  double getDerivative(const std::shared_ptr<std::vector<double>> &f,
                       const size_t &idx,
                       const Derivative &type);
  
};

class StructPropNew : public StructPropBase {

public:
  
  explicit StructPropNew(const VSStlsInput &in_);

private:

  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<std::shared_ptr<StlsCSRNew>> csr;
  // setup the csr vector
  std::vector<VSStlsInput> setupCSRInput(const VSStlsInput &in);
  void setupCSR(const VSStlsInput& in_);
  // 
  void doIterations();
  
  
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
