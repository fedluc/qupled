#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vsbase.hpp"
#include <limits>
#include <map>

class ThermoProp;
class StructProp;
class StlsCSRNew;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

class VSStls : public VSBase, public Stls {

public:

  // Constructor from initial data
  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in_);
  // Solve the scheme
  using VSBase::compute;

private:

  // Thermodynamic properties
  std::shared_ptr<ThermoProp> thermoProp;
  // Input parameters
  const VSInput &in() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
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

class StructProp : public StructPropBase {

public:

  explicit StructProp(const std::shared_ptr<const VSStlsInput> &in_);

private:

  // Structural properties
  std::shared_ptr<StlsCSRNew> csr;
  // Input parameters
  const VSStlsInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<IterationInput, VSStlsInput>(inPtr);
  }
};

// -----------------------------------------------------------------
// StlsCSRNew class
// -----------------------------------------------------------------

class StlsCSRNew : public CSRNew, public Stls {

public:

  // Constructor
  explicit StlsCSRNew(const std::shared_ptr<const VSStlsInput> &in_,
                      const bool isMaster_);
  // Compute static local field correction
  int compute() override;
  // // Getters
  const std::vector<double> &getSsf() const override { return Stls::getSsf(); }
  const std::vector<double> &getWvg() const override { return Stls::getWvg(); }
  const Vector2D &getLfc() const override { return Stls::getLfc(); }

private:

  // Flag marking wheter init was already called or not
  bool isInitialized;
  // Input parameters
  const VSInput &inVS() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
  }
  const Input &inRpa() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, Input>(inPtr);
  }
  // setup the csr vector
  void setupAuxiliaryStatePoints(const VSStlsInput &in);
  // Getters
  Vector2D &getLfc() override { return lfc; }
  void init() override;
  void initialGuessStls() override { Stls::initialGuess(); };
  void computeSsfStls() override { Stls::computeSsf(); };
  void computeLfc() override { CSRNew::computeLfc(); }
  void computeSsf() override { CSRNew::computeSsf();};
  void initialGuess() override { CSRNew::initialGuess(); };
  void computeLfcStls() override { Stls::computeLfc();};
};

#endif
