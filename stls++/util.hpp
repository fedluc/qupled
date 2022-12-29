#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <functional>
#include "numerics.hpp"

using namespace std;

namespace inpututil {

  // Types
  typedef const string cString;
  template<typename T> using cVector = const vector<T>;

  // Extract tokens from a string
  vector<string> tokenize(cString &str, const char separator);

  // Match a keyword from input to the underlying input data structures
  void matchKeyAndData(cVector<string> &keyword,
		       cString &input,
		       map<string,function<void(cString&, cString&)>> &funcArr);
  /// 
  void matchKeyAndData(cString &keyword,
		       cString &input,
		       map<string,function<void(cString&)>> &funcArr);
  
  // Compare numbers stored in strings
  template<typename T> bool isNegative(cString &str);
  ///
  template<typename T> bool isNotPositive(cString &str);
  ///
  template<typename T> bool isLarger(cString &str1, T num);
  
};

// Util functions to manipulate strings
namespace stringUtil{
  
  template<typename ... Args>
  string format(const string &format, Args ... args )
  {
    int size_s = snprintf( nullptr, 0, format.c_str(), args ... ) + 1;
    if( size_s <= 0 ) throw runtime_error( "Error during string formatting." );
    auto size = static_cast<size_t>( size_s );
    unique_ptr<char[]> buf( new char[ size ] );
    snprintf( buf.get(), size, format.c_str(), args ... );
    return string( buf.get(), buf.get() + size - 1 );
  }
  
}

// Util functions to manipulate vector
namespace vecUtil {

  // Element-wise sum between two vectors
  vector<double> sum(const vector<double> &v1,
		     const vector<double> &v2);
  // Element-wise difference between two vectors
  vector<double> diff(const vector<double> &v1,
		      const vector<double> &v2);
  // Element-wise multiplication of a vector and a scalar
  vector<double> mult(const vector<double> &v,
		      const double a);
  // Root mean square difference between two vectors
  double rms(const vector<double> &v1,
	     const vector<double> &v2,
	     const bool normalize);
}


namespace thermoUtil {

  // Class for internal energy calculation
  class InternalEnergy {
    
  private:

    // Coupling parameter
    const double rs = 0;
    // Integration limits
    const double yMin = 0;
    const double yMax = 0;
    // Integrator object
    const shared_ptr<Integrator1D> itg;
    // Static structure factor interpolator
    const shared_ptr<Interpolator> ssfi;
    // Integrand
    double integrand(const double y) const ;
    // Compute static structure factor
    double ssf(double x_) const;
    // Constant for unit conversion
    const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  
  public:

    // Constructor
    InternalEnergy(const double rs_,
		   const double yMin_,
		   const double yMax_,
		   const shared_ptr<Integrator1D> &itg_,
		   const shared_ptr<Interpolator> &ssfi_)
      : rs(rs_), yMin(yMin_), yMax(yMax_),
	itg(itg_), ssfi(ssfi_) {;};
    // Get result of integration 
    double get() const;
  
  };


  // Class for radial distribution function calculation
  class Rdf {
    
  private:

    // Spatial position
    const double r = 0;
    // Cutoff in the wave-vector grid
    const double cutoff = 0;
    // Integrator object
    const shared_ptr<Integrator1DFourier> itg;
    // Static structure factor interpolator
    const shared_ptr<Interpolator> ssfi;
    // Integrand
    double integrand(const double y) const ;
    // Compute static structure factor
    double ssf(double y_) const ;
  
  public:

    // Constructor
    Rdf(const double r_,
	const double cutoff_,
	const shared_ptr<Integrator1DFourier> &itg_,
	const shared_ptr<Interpolator> &ssfi_)
      : r(r_), cutoff(cutoff_), itg(itg_), ssfi(ssfi_) {
      assert(r>0);
      itg->setR(r);
    };
    // Get result of integration 
    double get() const;
  
  };
  
}
#endif
