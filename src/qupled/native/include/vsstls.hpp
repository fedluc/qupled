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
};

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

class StlsCSR : public CSR, public Stls {

public:

  // Constructor
  explicit StlsCSR(const std::shared_ptr<const VSStlsInput> &in_)
      : StlsCSR(in_, true) {}
  StlsCSR(const std::shared_ptr<const VSStlsInput> &in_, const bool isMaster_);
  // Solve the scheme
  int compute() override { return Stls::compute(); };

private:

  // Input parameters
  const VSInput &inVS() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
  }
  const Input &inRpa() const override {
    return *StlsUtil::dynamic_pointer_cast<Input, Input>(inPtr);
  }
  // setup the csr vector
  void setupWorkers(const VSStlsInput &in) {
    CSR::setupWorkers<StlsCSR, VSStlsInput>(in);
  }
  // Methods called by compute
  void init() override { CSR::init(); };
  void computeLfc() override { CSR::computeLfc(); };
  void computeSsf() override { CSR::computeSsf(); };
  void initialGuess() override { CSR::initialGuess(); };
  void updateSolution() override { CSR::updateSolution(); };
  double computeError() const override { return CSR::computeError(); }
  void initWorker() override { Stls::init(); };
  void computeLfcWorker() override { Stls::computeLfc(); };
  void computeSsfWorker() override { Stls::computeSsf(); };
  void initialGuessWorker() override { Stls::initialGuess(); };
  void updateSolutionWorker() override { Stls::updateSolution(); };
  double computeErrorWorker() const override { return Stls::computeError(); };
  Vector2D &getLfc() override { return lfc; };
  // Getters
  const std::vector<double> &getSsf() const override { return Stls::getSsf(); }
  const std::vector<double> &getWvg() const override { return Stls::getWvg(); }
  const Vector2D &getLfc() const override { return Stls::getLfc(); }
};

#endif
