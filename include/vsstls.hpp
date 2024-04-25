#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include <limits>
#include <map>
#include "vsbase.hpp"
#include "stls.hpp"

// Forward declarations
class VSStlsInput;

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

class StlsCSR : public CSR<std::vector<double>, Stls, VSStlsInput> {
    
private:

  using CSR = CSR<std::vector<double>, Stls, VSStlsInput>;
  // Helper methods to compute the derivatives
  double getDerivative(const std::shared_ptr<std::vector<double>>& f,
		       const size_t& idx,
		       const Derivative& type);
public:
  
  // Constructor
  explicit StlsCSR(const VSStlsInput& in_)
    : CSR(in_, Stls(in_, false, false)) { ; }
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();
  
};

class StructProp : public StructPropBase<StlsCSR, VSStlsInput> {

protected:

  using StructPropBase = StructPropBase<StlsCSR, VSStlsInput>;
  // Perform iterations to compute structural properties
  void doIterations();
  
public:

  // Constructor
  explicit StructProp(const VSStlsInput &in_) : StructPropBase(in_) { ; }
  
};


using ThermoProp = ThermoPropBase<StructProp, VSStlsInput>;

class VSStls : public VSBase<ThermoProp, Rpa, VSStlsInput> {

private:

  using VSBase = VSBase<ThermoProp, Rpa, VSStlsInput>;
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void updateSolution();
 
public:

  // Constructor from initial data
  explicit VSStls(const VSStlsInput &in_) : VSBase(in_) { ; }
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_,
	 const ThermoProp& thermoProp_) : VSBase(in_, thermoProp_) { ; }
  
};

#endif
