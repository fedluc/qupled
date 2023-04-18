#ifndef UTIL_HPP
#define UTIL_HPP

#include <fstream>
#include "numerics.hpp"

// Util functions to manipulate strings
namespace stringUtil {
  
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

  // Element-wise multiplication of two vectors
  vector<double> mult(const vector<double> &v1,
		      const vector<double> &v2);

  // Element-wise division of two vectors
  vector<double> div(const vector<double> &v1,
		     const vector<double> &v2);
  
  // Element-wise multiplication of a vector and a scalar
  vector<double> mult(const vector<double> &v,
		      const double a);

  // Root mean square difference between two vectors
  double rms(const vector<double> &v1,
	     const vector<double> &v2,
	     const bool normalize);

  // Class to represent 2D vectors
  class Vector2D {
  private:
    vector<double> v;
    size_t s1;
    size_t s2;
  public:
    Vector2D(const size_t s1_,
	     const size_t s2_)
      : v(s1_*s2_,0.0), s1(s1_), s2(s2_) {;};
    Vector2D()
      : Vector2D(0, 0) {;};
    size_t size() const;
    size_t size(const size_t i) const;
    void resize(const size_t s1_,
		const size_t s2_);
    double& operator()(const size_t i,
		       const size_t j);
    const double& operator()(const size_t i,
			     const size_t j) const;
    const double& operator()(const size_t i) const;
    vector<double>::iterator begin();
    vector<double>::iterator end();
    vector<double>::const_iterator begin() const;
    vector<double>::const_iterator end() const;
    void fill(const double &num);
    void fill(const size_t i,
	      const double &num);
    void fill(const size_t i,
	      const vector<double> &num);
    void sum(const Vector2D &v_);
    void diff(const Vector2D &v_);
    void mult(const Vector2D &v_);
    void mult(const double &num);
    void div(const Vector2D &v_);
  };
  
  // Class to represent 3D vectors
  class Vector3D {
  private:
    vector<double> v;
    size_t s1;
    size_t s2;
    size_t s3;
  public:
    Vector3D(const size_t s1_,
	     const size_t s2_,
	     const size_t s3_)
      : v(s1_*s2_*s3_,0.0), s1(s1_), s2(s2_), s3(s3_) {;};
    Vector3D()
      : Vector3D(0, 0, 0) {;};
    size_t size() const;
    size_t size(const size_t i) const;
    void resize(const size_t s1_,
		const size_t s2_,
		const size_t s3_);
    double& operator()(const size_t i,
		       const size_t j,
		       const size_t k);
    const double& operator()(const size_t i,
			     const size_t j,
			     const size_t k) const;
    const double& operator()(const size_t i,
			     const size_t j) const;
    const double& operator()(const size_t i) const;
    vector<double>::iterator begin();
    vector<double>::iterator end();
    vector<double>::const_iterator begin() const;
    vector<double>::const_iterator end() const;
    void fill(const double &num);
    void fill(const size_t i,
	      const size_t j,
	      const double &num);
    void fill(const size_t i,
	      const size_t j,
	      const vector<double> &num);
    void sum(const Vector3D &v_);
    void diff(const Vector3D &v_);
    void mult(const Vector3D &v_);
    void mult(const double &num);
    void div(const Vector3D &v_);
  };
  
    
}

// Util functions to compute thermodynamic properties
namespace thermoUtil {

  // Class for internal energy calculation
  class InternalEnergy {
    
  private:

    // Coupling parameter
    const double rs;
    // Integration limits
    const double yMin;
    const double yMax;
    // Integrator object
    Integrator1D &itg;
    // Static structure factor interpolator
    const Interpolator1D &ssfi;
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
		   const Interpolator1D &ssfi_,
		   Integrator1D &itg_) : rs(rs_), yMin(yMin_), yMax(yMax_),
					 itg(itg_), ssfi(ssfi_) {;};
    // Get result of integration 
    double get() const;
  
  };


  // Class for radial distribution function calculation
  class Rdf {
    
  private:

    // Spatial position
    const double r;
    // Cutoff in the wave-vector grid
    const double cutoff;
    // Integrator object
    Integrator1DFourier &itg;
    // Static structure factor interpolator
    const Interpolator1D &ssfi;
    // Integrand
    double integrand(const double y) const ;
    // Compute static structure factor
    double ssf(double y_) const ;
  
  public:

    // Constructor
    Rdf(const double r_,
	const double cutoff_,
	const Interpolator1D &ssfi_,
	Integrator1DFourier &itg_) : r(r_), cutoff(cutoff_), itg(itg_),
				     ssfi(ssfi_) { ; };
    // Get result of integration 
    double get() const;
  
  };

  double computeInternalEnergy(const vector<double> &wvg,
			       const vector<double> &ssf,
			       const double &coupling);
  
  vector<double> computeRdf(const vector<double> &r,
			    const vector<double> &wvg,
			    const vector<double> &ssf);
  
}

// Util functions to manipulate binary files
namespace binUtil {

  template<typename T>
  void writeNum(ofstream &file, const T &num) {
    file.write(reinterpret_cast<const char*>(&num), sizeof(num));
  };
  
  template<typename T>
  void writeDataToBinary(ofstream &file, const double &data) {
    writeNum<double>(file, data);
  };
  
  template<typename T>
  void writeDataToBinary(ofstream &file, const int &data) {
    writeNum<int>(file, data);
  };
  
  template<class T>
  void writeDataToBinary(ofstream &file, const T &data) {
    for (auto &el : data) { writeDataToBinary<decltype(el)>(file, el); }
  };
  
  template<typename T>
  void readNum(ifstream &file, T &num) {
    file.read((char*)&num, sizeof(T));
  };
  
  template<typename T>
  void readDataFromBinary(ifstream &file, double &data) {
    readNum<double>(file, data);
  };
  
  template<typename T>
  void readDataFromBinary(ifstream &file, int &data) {
    readNum<int>(file, data);
  };
  
  template<class T>
  void readDataFromBinary(ifstream &file, T &data) {
    for (auto &el : data) { readDataFromBinary<decltype(el)>(file, el);}
  };

}

#endif
