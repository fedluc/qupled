#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <functional>
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "util.hpp"

// -----------------------------------------------------------------
// C++ wrappers to GSL objects
// -----------------------------------------------------------------

namespace GslWrappers {

  // Wrapper to gsl_function
  template< typename T >
  class GslFunctionWrap : public gsl_function {

  private:

    const T& func;
    static double invoke(double x, void *params) {
      return static_cast<GslFunctionWrap*>(params)->func(x);
    }
  
  public:
  
    GslFunctionWrap(const T& func_) : func(func_) {
      function = &GslFunctionWrap::invoke;
      params = this;
    }
  
  };

  // Wrapper to handle errors in GSL functions
  template<typename Func, typename... Args>
  void callGSLFunction(Func&& gslFunction, Args&&... args);

  // Wrapper to handle allocation errors in GSL functions
  template<typename Ptr, typename Func, typename... Args>
  void callGSLAlloc(Ptr& ptr, Func&& gslFunction, Args&&... args);

}

// -----------------------------------------------------------------
// Classes to interpolate data
// -----------------------------------------------------------------

// Interpolator for 1D data
class Interpolator1D {

private:

  // Spline
  gsl_spline *spline;
  // Accelerator
  gsl_interp_accel *acc;
  // Cutoff (extrapolation for x > cutoff)
  double cutoff;
  // Size
  size_t n;
  // 
  // Setup interpolator
  void setup(const double &x,
	     const double &y,
	     const size_t n_);
  
public:

  // Constructor
  Interpolator1D(const std::vector<double> &x,
		 const std::vector<double> &y);
  Interpolator1D(const double &x,
		 const double &y,
		 const size_t n_);
  Interpolator1D();
  // Destructor
  ~Interpolator1D();
  // Reset
  void reset(const double &x,
	     const double &y,
	     const size_t n_);
  // Evaluate
  double eval(const double& x) const;
  
};

// Interpolator for 2D data
class Interpolator2D {

private:

  // Spline
  gsl_spline2d *spline;
  // Accelerator
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  // Size
  size_t nx;
  size_t ny;
  // Setup interpolator
  void setup(const double &x,
	     const double &y,
	     const double &z,
	     const int nx_,
	     const int ny_);
public:

  // Constructor
  Interpolator2D(const double &x,
		 const double &y,
		 const double &z,
		 const int nx_,
		 const int ny_);
  Interpolator2D(const Interpolator2D &it);
  Interpolator2D();
  // Destructor
  ~Interpolator2D();
  // Reset
  void reset(const double &x,
	     const double &y,
	     const double &z,
	     const int szx_,
	     const int szy_);
  // Evaluate
  double eval(const double& x,
	      const double& y) const;
};

// -----------------------------------------------------------------
// Classes to find roots of equations
// -----------------------------------------------------------------

class RootSolverBase {

protected:

  // Accuracy
  const double relErr;
  // Iterations
  const int maxIter;
  int iter;
  // Solver status
  int status;
  // Solution
  double sol;
  // Protected constructor
  RootSolverBase(const double& relErr_,
		 const int maxIter_) : relErr(relErr_),
				       maxIter(maxIter_),
				       iter(0),
				       status(GSL_CONTINUE) { ; };
  RootSolverBase() :  RootSolverBase(1.0e-10, 1000) { ; };
  
public:

  double getSolution() const { return sol; };
  
};

class BrentRootSolver : public RootSolverBase {

private:

  // Function to solve
  gsl_function *F;
  // Type of solver
  const gsl_root_fsolver_type *rst;
  // Solver
  gsl_root_fsolver *rs;
  
public:

  BrentRootSolver();
  ~BrentRootSolver();
  void solve(const std::function<double(double)>& func,
	     const std::vector<double>& guess);
};

class SecantSolver : public RootSolverBase {
  
public:

  SecantSolver(const double relErr_,
	       const int maxIter_) : RootSolverBase(relErr_, maxIter_) { ; };
  SecantSolver() { ; };
  void solve(const std::function<double(double)>& func,
	     const std::vector<double>& guess);
  
};

// -----------------------------------------------------------------
// Class to compute 1D integrals
// -----------------------------------------------------------------
  
class Integrator1D {

public:

  enum Type {
    DEFAULT,
    FOURIER,
    SINGULAR
  };
  
  class Param {
  public:
    using Limits = std::pair<double, double>;
    const Limits limits = Limits(numUtil::NaN, numUtil::NaN);
    const double fourierR = numUtil::NaN;
    Param(const Limits& limits_) : limits(limits_) { ; }
    Param(const double fourierR_) : fourierR(fourierR_) { ; }
  };
  
private:

  // Base class for all integrators derived from GSL
  class Base {
  protected:
    // Function to integrate
    gsl_function *F;
    // Integrator type
    const Type type;
    // Integration workspace limit
    const size_t limit;
    // Accuracy
    const double relErr;
    // Residual error
    double err;
    // Solution
    double sol;
  public:
    // Constructors
    Base(const Type& type_,
	 const size_t &limit_,
	 const double &relErr_) : type(type_),
				  limit(limit_),
				  relErr(relErr_) { ; }
    // Destructor
    virtual ~Base() = default;
    // Getters
    double getSolution() const { return sol; }
    double getAccuracy() const { return relErr; }
    Type getType() const { return type; }
    // Compute integral
    virtual void compute(const std::function<double(double)>& func,
			 const Param& param) = 0;
  
  };

  // CQUAD integrator from GSL
  class CQUAD : public Base {
  private:
    // Integration workspace
    gsl_integration_cquad_workspace *wsp;
    // Number of evaluations
    size_t nEvals;
  public:
    // Constructors
    CQUAD(const double &relErr_);
    CQUAD(const CQUAD& other) : Integrator1D::CQUAD(other.relErr) { ; }
    // Destructor
    ~CQUAD();
    // Compute integral
    void compute(const std::function<double(double)>& func,
		 const Param& param) override;
    
  };

  // QAWO integrator from GSL 
  class QAWO : public Base {
  private:
    // Integration workspace
    gsl_integration_workspace *wsp;
    gsl_integration_workspace *wspc;
    gsl_integration_qawo_table *qtab;
  public:
    // Constructors
    QAWO(const double& relErr_);
    QAWO(const QAWO& other) : Integrator1D::QAWO(other.relErr) { ; }
    // Destructor
    ~QAWO();
    // Compute integral
    void compute(const std::function<double(double)>& func,
		 const Param& param) override;

  };

  // QAGS integrator from GSL 
  class QAGS : public Base {
  private:
    // Integration workspace
    gsl_integration_workspace *wsp;
  public:
    // Constructors
    QAGS(const double& relErr_);
    QAGS(const QAGS& other) : Integrator1D::QAGS(other.relErr) { ; }
    // Destructor
    ~QAGS();
    // Compute integral
    void compute(const std::function<double(double)>& func,
		 const Param& param) override;
  };

  // Pointers to GSL integrals
  std::unique_ptr<Base> gslIntegrator;

public:

  // Constructors
  Integrator1D(const Type& type,
	       const double& relErr);
  Integrator1D(const double& relErr) : Integrator1D(Type::DEFAULT, relErr) { ; }
  Integrator1D(const Integrator1D& other) : Integrator1D(other.getType(),
							 other.getAccuracy()) { ; }
  // Compute integral
  void compute(const std::function<double(double)>& func,
	       const Param& param) const;
  // Getters
  double getSolution() const;
  double getAccuracy() const { return gslIntegrator->getAccuracy(); }
  Type getType() const { return gslIntegrator->getType(); }
  
};

// -----------------------------------------------------------------
// Class to compute 2D integrals
// -----------------------------------------------------------------

class Integrator2D {

private:

  // Typedef
  using Type = Integrator1D::Type;
  using Param = Integrator1D::Param;
  using Limits = Param::Limits;
  // Level 1 integrator (outermost integral)
  Integrator1D itg1;
  // Level 2 integrator
  Integrator1D itg2;
  // Temporary variable for level 2 integration
  double x;
  // Solution
  double sol;
  
public:

  // Constructors
  Integrator2D(const Type& type1,
	       const Type& type2,
	       const double& relErr) : itg1(type1, relErr),
				       itg2(type2, relErr) { ; }
  Integrator2D(const Type& type,
	       const double& relErr) : Integrator2D(type, type, relErr) { ; }
  Integrator2D(const double& relErr) : Integrator2D(Type::DEFAULT, relErr) { ; }
  // Compute integral
  void compute(const std::function<double(double)>& func1,
	       const std::function<double(double)>& func2,
	       const double& xMin,
	       const double& xMax,
	       const std::function<double(double)>& yMin,
	       const std::function<double(double)>& yMax,
	       const std::vector<double>& xGrid);
  void compute(const std::function<double(double)>& func1,
	       const std::function<double(double)>& func2,
	       const double& xMin,
	       const double& xMax,
	       const double& yMin,
	       const double& yMax,
	       const std::vector<double>& xGrid);
  // Getters
  double getX() const { return x; };
  double getSolution() const { return sol; };
};


#endif
