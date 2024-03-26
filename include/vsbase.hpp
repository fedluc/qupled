#ifndef VSBASE_HPP
#define VSBASE_HPP

#include <limits>
#include <map>

// Forward declarations
class VSStlsInput;
class SecantSolver;
class Interpolator1D;

// -----------------------------------------------------------------
// VSBase class
// -----------------------------------------------------------------

template<typename ThermoProp, typename Scheme, typename Input>
class VSBase : public Scheme {

protected: 

  // Input data
  Input in;
  // Thermodynamic properties
  ThermoProp thermoProp;
  // Free parameter
  double alpha;
  // Output verbosity
  const bool verbose;
  
  // Compute free parameter
  virtual double computeAlpha() = 0;
  
  // Iterations to solve the vs scheme
  void doIterations() {
    auto func = [this](const double& alphaTmp)->double{return alphaDifference(alphaTmp);};
    SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
    rsol.solve(func, in.getAlphaGuess());
    // if (!rsol.success()) {
    //   parallelUtil::MPI::throwError("VSStls: the root solver did not converge "
		// 		    "to the desired accuracy.");
    // }
    alpha = rsol.getSolution();
    if (verbose) { std::cout << "Free parameter = " << alpha << std::endl; }
    updateSolution();
  }

  // Object function used in the secant solver
  double alphaDifference(const double& alphaTmp) {
    alpha = alphaTmp;
    thermoProp.setAlpha(alpha);
    const double alphaTheoretical = computeAlpha();
    return alpha - alphaTheoretical;
  }
  
  // Update structural output solution
  virtual void updateSolution() = 0;

 
public:

  // Constructor from initial data
  VSBase(const Input &in_) : Scheme(in_),
			     in(in_),
			     thermoProp(in_),
			     verbose(true) { ; }
  // Constructor for recursive calculations
  VSBase(const Input &in_,
	 const ThermoProp& thermoProp_) : Scheme(in_),
					  in(in_),
					  thermoProp(in_, thermoProp_),
					  verbose(false) { ; }

  // Destructor
  virtual ~VSBase() = default;
  
  // Compute vs-stls scheme
  int compute() {
    try {
      Scheme::init();
      if (verbose) std::cout << "Free parameter calculation ..." << std::endl;
      doIterations();
      if (verbose) std::cout << "Done" << std::endl;
      return 0;
    }
    catch (const std::runtime_error& err) {
      std::cerr << err.what() << std::endl;
      return 1;
    }
  }
  
  // Getters
  const ThermoProp& getThermoProp() const {
    return thermoProp;
  }
  
  std::vector<std::vector<double>> getFreeEnergyIntegrand() const {
    return thermoProp.getFreeEnergyIntegrand();
  }
  
  std::vector<double> getFreeEnergyGrid() const {
    return thermoProp.getFreeEnergyGrid();
  }
  
};

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

template<typename StructProp, typename Input>
class ThermoPropBase {

protected:

  using SIdx = typename StructProp::Idx;
  enum Idx {
    THETA_DOWN,
    THETA,
    THETA_UP
  };
  // Map between struct and thermo indexes
  static constexpr int NPOINTS = 3;
  // Structural properties
  StructProp structProp;
  // Grid for thermodyamic integration
  std::vector<double> rsGrid;
  // Free energy integrand for NPOINTS state points
  std::vector<std::vector<double>> fxcIntegrand;
  // Flags marking particular state points
  bool isZeroCoupling;
  bool isZeroDegeneracy;
  // Compute the free energy
  double computeFreeEnergy(const SIdx iStruct,
			   const bool normalize) const {
    Idx iThermo;
    switch (iStruct) {
    case SIdx::RS_DOWN_THETA_DOWN: case SIdx::RS_THETA_DOWN: case SIdx::RS_UP_THETA_DOWN:
      iThermo = THETA_DOWN; break;
    case SIdx::RS_DOWN_THETA: case SIdx::RS_THETA: case SIdx::RS_UP_THETA:
      iThermo = THETA; break;
    case SIdx::RS_DOWN_THETA_UP: case SIdx::RS_THETA_UP: case SIdx::RS_UP_THETA_UP:
      iThermo = THETA_UP; break;
    default:
      assert(false);
      iThermo = THETA;
      break;
    }
    const std::vector<double>& rs = structProp.getCouplingParameters();
    return thermoUtil::computeFreeEnergy(rsGrid, fxcIntegrand[iThermo], rs[iStruct], normalize);
  }
  
public:

  // Constructors
  ThermoPropBase(const Input &in) : structProp(in) {
    const double& rs = in.getCoupling();
    const double& drs = in.getCouplingResolution();
    // Check if we are solving for particular state points
    isZeroCoupling = (rs == 0.0);
    isZeroDegeneracy = (in.getDegeneracy() == 0.0);
    // Build integration grid
    rsGrid.push_back(0.0);
    const double rsMax = rs + drs;
    while(!numUtil::equalTol(rsGrid.back(), rsMax)){
      rsGrid. push_back(rsGrid.back() + drs);
    }
    // Initialize the free energy integrand
    fxcIntegrand.resize(NPOINTS);
    const size_t nrs = rsGrid.size();
    for (auto& f : fxcIntegrand) {
      f.resize(nrs);
      vecUtil::fill(f, numUtil::Inf);
    }
    // Fill the free energy integrand if passed in input
    const auto& fxciData = in.getFreeEnergyIntegrand();
    if (!fxciData.grid.empty()) {
      for (const auto& theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
	const Interpolator1D itp(fxciData.grid, fxciData.integrand[theta]);
	const double rsMaxi = fxciData.grid.back();
	for (size_t i = 0; i < nrs; ++i) {
	  const double& rs = rsGrid[i];
	  if (rs <= rsMaxi) { fxcIntegrand[theta][i] = itp.eval(rs); }
	}
      }
    }
  }
  
  ThermoPropBase(const Input &in,
		 const ThermoPropBase &other) : ThermoPropBase(in) {
    assert(other.rsGrid[1] - other.rsGrid[0] == rsGrid[1] - rsGrid[0]);
    const size_t nrs = rsGrid.size();
    const double rsMax = other.rsGrid.back();
    for (const auto& theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
      const auto& fxciBegin = fxcIntegrand[theta].begin();
      const auto& fxciEnd = fxcIntegrand[theta].end();
      const auto& it = std::find(fxciBegin, fxciEnd, numUtil::Inf);
      size_t i =  std::distance(fxciBegin, it);
      while (i < nrs && rsGrid[i] < rsMax) {
	fxcIntegrand[theta][i] = other.fxcIntegrand[theta][i];
	++i;
      }
    }
  }
  
  // Set the value of the free parameter in the structural properties
  void setAlpha(const double& alpha) {
    structProp.setAlpha(alpha);
  }
  
  // Compute the thermodynamic properties
  template<typename Scheme>
  void compute(const Input &in) {
    // Recursive calls to solve the VS-STLS scheme for all state points
    // with coupling parameter smaller than rs
    const double nrs = rsGrid.size();
    Input inTmp = in;
    std::vector<double> fxciTmp(StructProp::NPOINTS);
    for (size_t i = 0; i < nrs; ++i) {
      const double& rs = rsGrid[i];
      if (numUtil::equalTol(rs, in.getCoupling())) {
	structProp.compute();
	fxciTmp = structProp.getFreeEnergyIntegrand();
      }
      else if (rs < in.getCoupling()) {
	if (rs == 0.0 || fxcIntegrand[THETA][i] != numUtil::Inf) { continue; }
	printf("Free energy integrand calculation, solving VS scheme for rs = %.5f:\n", rs);
	inTmp.setCoupling(rs);
	Scheme schemeTmp(inTmp, *this);
	schemeTmp.compute();
	fxciTmp = schemeTmp.getThermoProp().structProp.getFreeEnergyIntegrand();
	printf("Done\n");
	printf("---------------------------------------------------------------------------\n");
      }
      else {
	break;
      }
      fxcIntegrand[THETA_DOWN][i-1] = fxciTmp[SIdx::RS_DOWN_THETA_DOWN];
      fxcIntegrand[THETA_DOWN][i]   = fxciTmp[SIdx::RS_THETA_DOWN];
      fxcIntegrand[THETA_DOWN][i+1] = fxciTmp[SIdx::RS_UP_THETA_DOWN];
      fxcIntegrand[THETA][i-1]      = fxciTmp[SIdx::RS_DOWN_THETA];
      fxcIntegrand[THETA][i]        = fxciTmp[SIdx::RS_THETA];
      fxcIntegrand[THETA][i+1]      = fxciTmp[SIdx::RS_UP_THETA];
      fxcIntegrand[THETA_UP][i-1]   = fxciTmp[SIdx::RS_DOWN_THETA_UP];
      fxcIntegrand[THETA_UP][i]     = fxciTmp[SIdx::RS_THETA_UP];
      fxcIntegrand[THETA_UP][i+1]   = fxciTmp[SIdx::RS_UP_THETA_UP];
    }
  }

  // Get structural properties
  template<typename CSR>
  const CSR& getStructProp() {
    if (!structProp.isComputed()) { structProp.compute(); }
    if (isZeroCoupling && isZeroDegeneracy) {
      return structProp.getCsr(SIdx::RS_DOWN_THETA_DOWN); 
    }
    if (!isZeroCoupling && isZeroDegeneracy) {
      return structProp.getCsr(SIdx::RS_THETA_DOWN); 
    }
    if (isZeroCoupling && !isZeroDegeneracy) {
      return structProp.getCsr(SIdx::RS_DOWN_THETA); 
    }
    return structProp.getCsr(SIdx::RS_THETA); 
  }
  
  // Get free energy and free energy derivatives
  std::vector<double> getFreeEnergyData() const {
    const std::vector<double> rsVec = structProp.getCouplingParameters();
    const std::vector<double> thetaVec = structProp.getDegeneracyParameters();
    // Free energy
    const double fxc = computeFreeEnergy(SIdx::RS_THETA, true);
    // Free energy derivatives with respect to the coupling parameter
    double fxcr;
    double fxcrr;
    {
      const double rs = rsVec[SIdx::RS_THETA];
      const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
      const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA, false);
      const double f1 = computeFreeEnergy(SIdx::RS_THETA, false);
      const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA, false);
      fxcr = (f0 - f2) / (2.0 * drs * rs) - 2.0 * fxc;
      fxcrr = (f0 - 2.0 * f1 + f2) / (drs * drs) - 2.0 * fxc - 4.0 * fxcr;
    }
    // Free energy derivatives with respect to the degeneracy parameter
    double fxct;
    double fxctt;
    {
      const double theta = thetaVec[SIdx::RS_THETA];
      const double theta2 = theta * theta;
      const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
      const double f0 = computeFreeEnergy(SIdx::RS_THETA_UP, true);
      const double f1 = computeFreeEnergy(SIdx::RS_THETA_DOWN, true);
      fxct = theta * (f0 - f1) / (2.0 * dt);
      fxctt = theta2 * (f0 - 2.0 * fxc + f1) / (dt * dt);
    }
    // Free energy mixed derivatives
    double fxcrt;
    {
      const double t_rs = thetaVec[SIdx::RS_THETA] / rsVec[SIdx::RS_THETA];
      const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
      const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
      const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA_UP, false);
      const double f1 = computeFreeEnergy(SIdx::RS_UP_THETA_DOWN, false);
      const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA_UP, false);
      const double f3 = computeFreeEnergy(SIdx::RS_DOWN_THETA_DOWN, false);
      fxcrt = t_rs * (f0 - f1 - f2 + f3) / (4.0 * drs * dt) - 2.0 * fxct;
    }
    return std::vector<double>({fxc, fxcr, fxcrr, fxct, fxctt, fxcrt});
  }
  
  // Get internal energy and internal energy derivatives
  std::vector<double> getInternalEnergyData() const {
    // Internal energy
    const std::vector<double> uVec = structProp.getInternalEnergy();
    const double u = uVec[SIdx::RS_THETA];
    // Internal energy derivative with respect to the coupling parameter
    double ur;
    {
      const std::vector<double> rs = structProp.getCouplingParameters(); 
      const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
      const std::vector<double> rsu = structProp.getFreeEnergyIntegrand();
      const double& u0 = rsu[SIdx::RS_UP_THETA];
      const double& u1 = rsu[SIdx::RS_DOWN_THETA];
      ur = (u0 - u1) / (2.0 * drs) - u;
    }
    // Internal energy derivative with respect to the degeneracy parameter
    double ut;
    {
      const std::vector<double> theta = structProp.getDegeneracyParameters();
      const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
      const double u0 = uVec[SIdx::RS_THETA_UP];
      const double u1 = uVec[SIdx::RS_THETA_DOWN];
      ut = theta[SIdx::RS_THETA] * (u0 - u1) / (2.0 * dt);
    }
    return std::vector<double>({u, ur, ut});
  }
  
  // Get free energy integrand
  const std::vector<std::vector<double>>& getFreeEnergyIntegrand() const {
    return fxcIntegrand;
  }
  
  // Get free energy grid
  const std::vector<double>& getFreeEnergyGrid() const  {
    return rsGrid;
  }
  
};

// -----------------------------------------------------------------
// StructPropBase class
// -----------------------------------------------------------------

template<typename CSR, typename Input>
class StructPropBase {

public:
  
  enum Idx {
    RS_DOWN_THETA_DOWN,
    RS_THETA_DOWN,
    RS_UP_THETA_DOWN,
    RS_DOWN_THETA,
    RS_THETA,
    RS_UP_THETA,
    RS_DOWN_THETA_UP,
    RS_THETA_UP,
    RS_UP_THETA_UP,
  };
  static constexpr int NRS = 3;
  static constexpr int NTHETA = 3;
  static constexpr int NPOINTS = NRS * NTHETA;
  
protected:

  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<CSR> csr;
  // Flag marking whether the initialization for the stls data is done
  bool csrIsInitialized;
  // Flag marking whether the structural properties were computed
  bool computed;
  // Vector used as output parameter in the getters functions
  mutable std::vector<double> outVector;

  // Perform iterations to compute structural properties
  virtual void doIterations() = 0;

  // Generic getter function to return vector data
  const std::vector<double>& getBase(std::function<double(const CSR&)> f) const {
    for (size_t i = 0; i < NPOINTS; ++i) {
      outVector[i] = f(csr[i]);
    }
    return outVector; 
  }

  
public:

  // Constructor
  StructPropBase(const Input &in_) : csrIsInitialized(false),
				     computed(false),
				     outVector(NPOINTS) {
    Input inTmp = in_;
    const double& drs = inTmp.getCouplingResolution();
    const double& dTheta = inTmp.getDegeneracyResolution();
    // If there is a risk of having negative state parameters, shift the
    // parameters so that rs - drs = 0 and/or theta - dtheta = 0
    if (inTmp.getCoupling() < drs) { inTmp.setCoupling(drs); }
    if (inTmp.getDegeneracy() < dTheta) { inTmp.setDegeneracy(dTheta); }
    double rs = inTmp.getCoupling();
    double theta = inTmp.getDegeneracy();
    // Setup objects
    for (const double& thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
      for (const double& rsTmp : {rs - drs, rs, rs + drs}){
	inTmp.setDegeneracy(thetaTmp);
	inTmp.setCoupling(rsTmp);
	csr.push_back(CSR(inTmp));
      }
    }
    assert(csr.size() == NPOINTS);
    // Setup derivative dependency in the CSR objects  
    for (size_t i = 0; i < csr.size(); ++i) {
      switch (i) {
      case RS_DOWN_THETA_DOWN: case RS_DOWN_THETA: case RS_DOWN_THETA_UP:
	csr[i].setDrsData(csr[i + 1], csr[i + 2],
			  CSR::Derivative::FORWARD); break;
      case RS_THETA_DOWN: case RS_THETA: case RS_THETA_UP:
	csr[i].setDrsData(csr[i + 1], csr[i - 1],
			  CSR::Derivative::CENTERED); break;
      case RS_UP_THETA_DOWN: case RS_UP_THETA: case RS_UP_THETA_UP:
	csr[i].setDrsData(csr[i - 1], csr[i - 2],
			  CSR::Derivative::BACKWARD); break;
      }
    }
    for (size_t i = 0; i < csr.size(); ++i) {
      switch (i) {
      case RS_DOWN_THETA_DOWN: case RS_THETA_DOWN: case RS_UP_THETA_DOWN:
	csr[i].setDThetaData(csr[i + NRS], csr[i + 2 * NRS],
			     CSR::Derivative::FORWARD); break;
      case RS_DOWN_THETA: case RS_THETA: case RS_UP_THETA:
	csr[i].setDThetaData(csr[i + NRS], csr[i - NRS],
			     CSR::Derivative::CENTERED); break;
      case RS_DOWN_THETA_UP: case RS_THETA_UP: case RS_UP_THETA_UP:
	csr[i].setDThetaData(csr[i - NRS], csr[i - 2 * NRS],
			     CSR::Derivative::BACKWARD); break;
      }
    }
  }
  
  // Compute structural properties
  int compute() {
    try {
      if (!csrIsInitialized) {
	for (auto& c : csr) { c.init(); }
	csrIsInitialized = true;
      }
      doIterations();
      computed = true;
      return 0;
    }
    catch (const std::runtime_error& err) {
      std::cerr << err.what() << std::endl;
      return 1;
    }
  }
  
  // Set free parameter
  void setAlpha(const double& alpha) {
    for (auto& c : csr) { c.setAlpha(alpha); }
  }
  
  // Get coupling parameters for all the state points
  std::vector<double> getCouplingParameters() const {
    return getBase([&](const CSR& c){ return c.getCoupling();});
  }
  // Get degeneracy parameters for all the state points
  std::vector<double> getDegeneracyParameters() const {
    return getBase([&](const CSR& c){return c.getDegeneracy();});
  }
  
  // Get internal energy for all the state points
  std::vector<double> getInternalEnergy() const {
    return getBase([&](const CSR& c){return c.getInternalEnergy();});
  }
  
  // Get free energy integrand for all the state points
  std::vector<double> getFreeEnergyIntegrand() const{
    return getBase([&](const CSR& c){return c.getFreeEnergyIntegrand();});
  }
  
  // Get structural properties for output
  const CSR& getCsr(const Idx& idx) const { return csr[idx]; }
  
  // Boolean marking whether the structural properties where computed or not
  bool isComputed() const { return computed; }
  
};

// -----------------------------------------------------------------
// CSR base class
// -----------------------------------------------------------------

template <typename T, typename Scheme, typename Input>
class CSR : public Scheme {

public:

  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED,
		    FORWARD,
		    BACKWARD };
  
protected:
  
  // Default value of alpha
  static constexpr double DEFAULT_ALPHA = numUtil::Inf;
  // Input data
  const Input in;
  // local field correction (static or dynamic)
  T lfc;
  // Free parameter
  double alpha;
  // Pointer to the local field correction with rs = rs + drs
  T* lfcRsUp;
  // Pointer to the local field correction with rs = rs - drs
  T* lfcRsDown;
  // Pointer to the local field correction with theta = theta + dtheta
  T* lfcThetaUp;
  // Pointer to the local field correction with theta = theta - dtheta
  T* lfcThetaDown;
  // Numerical scheme used to compute the coupling parameter derivative
  Derivative dTypeRs;
  // Numerical scheme used to compute the degeneracy parametere derivative
  Derivative dTypeTheta;
    
  // Helper methods to compute the derivatives
  double getDerivative(const double& f0,
		       const double& f1,
		       const double& f2,
		       const Derivative& type) {
    switch(type) {
    case BACKWARD:
      return 3.0 * f0 - 4.0 * f1 + f2; break;
    case CENTERED:
      return f1 - f2; break;
    case FORWARD:
      return -getDerivative(f0, f1, f2, BACKWARD); break;
    default:
      assert(false);
      return -1;
      break;
    }
  }
  
public:

  // Constructor
  CSR(const Input& in_,
      const Scheme& scheme) : Scheme(scheme),
			      in(in_),
			      alpha(DEFAULT_ALPHA),
			      lfcRsUp(nullptr),
			      lfcRsDown(nullptr),
			      lfcThetaUp(nullptr),
			      lfcThetaDown(nullptr),
			      dTypeRs(CENTERED),
			      dTypeTheta(CENTERED) { ; }

  // Set the data to compute the coupling parameter derivative
  void setDrsData(CSR<T, Scheme, Input> &csrRsUp,
		  CSR<T, Scheme, Input> &csrRsDown,
		  const Derivative &dTypeRs) {
    this->lfcRsUp = &csrRsUp.lfc;
    this->lfcRsDown = &csrRsDown.lfc;
    this->dTypeRs = dTypeRs;
  }
  
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(CSR<T, Scheme, Input> &csrThetaUp,
		     CSR<T, Scheme, Input> &csrThetaDown,
		     const Derivative &dTypeTheta) {
    this->lfcThetaUp = &csrThetaUp.lfc;
    this->lfcThetaDown = &csrThetaDown.lfc;
    this->dTypeTheta = dTypeTheta;
  }

  // Initialize the scheme
  void init() {
    Scheme::init();
  }
  
  // Set the free parameter
  void setAlpha(const double& alpha) {
    this->alpha = alpha;
  }

  // Get the coupling parameter
  double getCoupling() const {
    return in.getCoupling();
  }

  // Get the degenarcy parameter
  double getDegeneracy() const {
    return in.getDegeneracy();
  }
  
  // Compute the internal energy
  double getInternalEnergy() const {
    const double rs = in.getCoupling();
    return thermoUtil::computeInternalEnergy(Scheme::wvg, Scheme::ssf, rs);
  }
  
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const {
    return thermoUtil::computeInternalEnergy(Scheme::wvg, Scheme::ssf, 1.0);
  }
  

};

#endif
