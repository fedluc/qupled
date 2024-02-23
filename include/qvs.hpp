#ifndef QVS_HPP
#define QVS_HPP

#include <limits>
#include <map>
#include <functional>
#include "vsstls.hpp"
#include "qstls.hpp"


// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class qVSStls : public VSStls, public Qssf {

protected:

    void doIterations();
    qThermoProp qthermoProp;
    double computeAlpha();
    double alphaDifference(const double& alphaTmp);
    

public:
    void updateSolution();
    int compute();
  
};


// -----------------------------------------------------------------
// Class to handle simultaneous state point calculations
// -----------------------------------------------------------------


template<typename T = std::vector<double>>
class qStlsCSR : public CSR<T>, public Qstls {

    friend class qStructProp;

protected:

  void computeAdrStls();
  void computeAdr();
  // Set the free parameter
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  // Compute Q
  double getQ() const;
  

public:
  // Constructor
  qStlsCSR(const QstlsInput& in_) : Qstls(in_),
				    in(in_),
				    alpha(DEFAULT_ALPHA),
				    stlsAdrRsUp(nullptr),
				    stlsAdrRsDown(nullptr),
				    stlsAdrThetaUp(nullptr),
				    stlsAdrThetaDown(nullptr),
				    dTypeRs(CENTERED),
				    dTypeTheta(CENTERED) { ; }
  
  
};


// -----------------------------------------------------------------
// Class to handle the state point derivatives
// -----------------------------------------------------------------

class qStructProp : public StructProp {

protected:
    void doIterations();
    // Flag marking whether the initialization for the stls data is done
    bool stlsAdrIsInitialized;
    // Vector containing NPOINTS state points to be solved simultaneously
    std::vector<qStlsCSR<std::vector<double>>> Adr;
    // Get Q
    
    std::vector<double> getCouplingParameters() const;
    std::vector<double> getDegeneracyParameters() const;
    std::vector<double> getInternalEnergy() const;
    std::vector<double> getFreeEnergyIntegrand() const;
    // Generic getter function to return vector data
    const std::vector<double>& getBase(std::function<double(const qStlsCSR<std::vector<double>>&)> f) const;

    template<typename T = std::vector<double>>
    const qStlsCSR<T>& getAdrStls(const Idx& idx) const { return Adr[idx]; }

    // Flag marking whether the structural properties were computed
    bool computed;

public:

    int compute();
    void setAlpha(const double& alpha);
    std::vector<double> getQ() const;
    // Boolean marking whether the structural properties where computed or not
    bool isComputed() const { return computed; }

};


// -----------------------------------------------------------------
// Class to handle the thermodynamic properties and derivatives
// -----------------------------------------------------------------

class qThermoProp : public ThermoProp {

protected:

    // Compute the free energy
    double computeQ(const SIdx iStruct,
			   const bool normalize) const ;
    qStructProp qstructProp;
    void setAlpha(const double& alpha);

    

public:

    void compute(const VSStlsInput& in);
    // Get internal energy and internal energy derivatives
    std::vector<double> getQData() const;
    
    template<typename T = std::vector<double>>
    const qStlsCSR<T>& getStructProp();
};

class Q : public AdrFixed {
    // Get integration result
    void get(std::vector<double> &wvg,
	   vecUtil::Vector3D &res) const;
    double integranddenom(const double q,
		    const double l) const;
    double integrandnum(const double t,
		    const double y,
		    const double l) const;

};


#endif