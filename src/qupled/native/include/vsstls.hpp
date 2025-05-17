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
      : VSStls(std::make_shared<const VSStlsInput>(in_)){};
  // Solve the scheme
  using VSBase::compute;

private:

  // Thermodynamic properties
  std::shared_ptr<ThermoProp> thermoProp;
  // Input parameters
  const VSInput &in() const override {
    return *StlsUtil::castInput<Input, VSInput>(inPtr);
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
  explicit StlsCSR(const std::shared_ptr<const VSStlsInput> &in_)
      : CSR(),
        Stls(in_, false) {}
  // Compute static local field correction
  void computeLfcStls();
  void computeLfc() override;
  // Publicly esposed private stls methods
  void init() override { Stls::init(); }
  void initialGuess() override { Stls::initialGuess(); }
  void computeSsf() override { Stls::computeSsf(); }
  double computeError() override { return Stls::computeError(); }
  void updateSolution() override { Stls::updateSolution(); }
  // Getters
  const std::vector<double> &getSsf() const override { return Stls::getSsf(); }
  const std::vector<double> &getWvg() const override { return Stls::getWvg(); }
  const Vector2D &getLfc() const override { return Stls::getLfc(); }

private:

  // Input parameters
  const VSInput &inVS() const override {
    return *StlsUtil::castInput<Input, VSInput>(inPtr);
  }
  const Input &inRpa() const override {
    return *StlsUtil::castInput<Input, Input>(inPtr);
  }
};

#endif
