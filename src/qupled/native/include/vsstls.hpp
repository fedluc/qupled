#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "esa.hpp"
#include "input.hpp"
#include "stls.hpp"
#include "vsbase.hpp"
#include <limits>
#include <map>

class ThermoProp;
class StructProp;
class StlsCSR;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

class VSStls : public VSBase, public Stls {

public:

  // Constructor from initial data
  explicit VSStls(const VSStlsInput &in_);
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_, const ThermoProp &thermoProp_);
  // Solve the scheme
  using VSBase::compute;
  // Getters
  const VSStlsInput &getInput() const { return in; }

private:

  // Input
  VSStlsInput in;
  // Thermodynamic properties
  std::shared_ptr<ThermoProp> thermoProp;
  // Initialize
  void initScheme();
  void initFreeEnergyIntegrand();
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
  // Print info
  void print(const std::string &msg) { VSBase::print(msg); }
  void println(const std::string &msg) { VSBase::println(msg); }
};

// -----------------------------------------------------------------
// ThermoProp class
// -----------------------------------------------------------------

class ThermoProp : public ThermoPropBase {

public:

  // Constructor
  explicit ThermoProp(const VSStlsInput &in_);

private:

  // Structural properties
  std::shared_ptr<StructProp> structProp;
};

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

class StructProp : public Logger, public StructPropBase {

public:

  explicit StructProp(const VSStlsInput &in_);

private:

  // Input
  const VSStlsInput in;
  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<std::shared_ptr<StlsCSR>> csr;
  // setup the csr vector
  std::vector<VSStlsInput> setupCSRInput();
  void setupCSR();
  //
  void doIterations();
};

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

class StlsCSR : public CSR, public ESA {

public:

  // Constructor
  explicit StlsCSR(const VSStlsInput &in_)
      : CSR(in_, in_),
        ESA(in_, false),
        in(in_) {}

  // Compute static local field correction
  void computeSlfc() { ESA::computeSlfc(); };

  // Publicly esposed private stls methods
  void init() { ESA::init(); }
  void initialGuess() {}
  void computeSsf() { ESA::computeSsf(); }
  double computeError() { return 0; }
  void updateSolution() {}

  // Getters
  const std::vector<double> &getSsf() const { return ESA::getSsf(); }
  const std::vector<double> &getSlfc() const { return ESA::getSlfc(); }
  const std::vector<double> &getWvg() const { return ESA::getWvg(); }

private:

  // Input parameters
  VSStlsInput in;
};

#endif
