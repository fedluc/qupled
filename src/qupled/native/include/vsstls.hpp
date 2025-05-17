#ifndef VSSTLS_HPP
#define VSSTLS_HPP

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
  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in_);
  explicit VSStls(const VSStlsInput &in_)
      : VSStls(std::make_shared<const VSStlsInput>(in_)) {};
  // Solve the scheme
  using VSBase::compute;

private:

  // Thermodynamic properties
  std::shared_ptr<ThermoProp> thermoProp;
  // Input parameters
  const VSInput &in() const override {
    return dynamic_cast<const VSInput &>(*inPtr);
  }
  // Initialize
  void init() override;
  // Compute free parameter
  double computeAlpha() override;
  // Iterations to solve the vs-stls scheme
  void updateSolution() override;
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
  explicit ThermoProp(const std::shared_ptr<const VSStlsInput> &in_);

private:

  // Structural properties
  std::shared_ptr<StructProp> structProp;
};

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

class StructProp : public Logger, public StructPropBase {

public:

  explicit StructProp(const std::shared_ptr<const VSStlsInput> &in_);

private:

  // Input
  const std::shared_ptr<const VSStlsInput> in;
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

class StlsCSR : public CSR, public Stls {

public:

  // Constructor
  explicit StlsCSR(const VSStlsInput &in_)
      : CSR(in_, in_),
        Stls(in_, false),
        in(in_) {}

  // Compute static local field correction
  void computeLfcStls();
  void computeLfc();

  // Publicly esposed private stls methods
  void init() { Stls::init(); }
  void initialGuess() { Stls::initialGuess(); }
  void computeSsf() { Stls::computeSsf(); }
  double computeError() { return Stls::computeError(); }
  void updateSolution() { Stls::updateSolution(); }

  // Getters
  const std::vector<double> &getSsf() const { return Stls::getSsf(); }
  const std::vector<double> &getWvg() const { return Stls::getWvg(); }
  const Vector2D &getLfc() const { return Stls::getLfc(); }

private:

  // Input parameters
  VSStlsInput in;
};

#endif
