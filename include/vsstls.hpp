#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vsbase.hpp"
#include <limits>
#include <map>

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

class StlsCSR : public Stls, public CSRNew {

  // Typedef
  using Lfc = std::vector<double>;
  using LfcPtr = std::shared_ptr<Lfc>;
  
public:

  // Data for the local field correction with modified state point
  struct DerivativeData {
    Derivative type;
    LfcPtr up;
    LfcPtr down;
  };
  
  // Constructor
  explicit StlsCSR(const VSStlsInput &in_)
    : Stls(in_.toStlsInput(), false, false),
      in(in_),
      lfc(std::make_shared<Lfc>()) {}
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();

  // Publicly exposed private stls methods
  void init() { Stls::init(); }
  void initialGuess() { Stls::initialGuess(); }
  void computeSsf() { Stls::computeSsf(); }
  double computeError() { return Stls::computeError(); }
  void updateSolution() { Stls::updateSolution(); }
  std::vector<double> getWvg() const { return Stls::getWvg(); }
  std::vector<double> getSsf() const { return Stls::getSsf(); }
  VSStlsInput const getInput() const { return in; }
  double getCoupling() const { return in.getCoupling(); }
  double getDegeneracy() const { return in.getDegeneracy(); }

private:

  // Input
  VSStlsInput in;
  // static local field correction
  LfcPtr lfc;
  // Data for the local field correction with modified coupling paramter
  DerivativeData lfcRs;
  // Data for the local field correction with modified degeneracy parameter
  DerivativeData lfcTheta;
  
  // Helper methods to compute the derivatives
  double getDerivative(const std::shared_ptr<std::vector<double>> &f,
                       const size_t &idx,
                       const Derivative &type);
};

class StructProp : public StructPropBase<StlsCSR, VSStlsInput> {

public:

  // Constructor
  explicit StructProp(const VSStlsInput &in_)
      : StructPropBase(in_) {}

protected:

  using StructPropBase = StructPropBase<StlsCSR, VSStlsInput>;
  // Perform iterations to compute structural properties
  void doIterations();
};

using ThermoProp = ThermoPropBase<StructProp, VSStlsInput>;

class VSStls : public VSBase<ThermoProp, Rpa, VSStlsInput> {

public:

  // Constructor from initial data
  explicit VSStls(const VSStlsInput &in_)
      : VSBase(in_) {}
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_, const ThermoProp &thermoProp_)
      : VSBase(in_, thermoProp_) {}

private:

  using VSBase = VSBase<ThermoProp, Rpa, VSStlsInput>;
  // Compute free parameter
  double computeAlpha() override;
  // Iterations to solve the vs-stls scheme
  void updateSolution() override;
  // Fill the free energy integrand
  void fillFreeEnergyIntegrand() override;
};

#endif
