#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
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
  template<typename T> bool isLarger(cString &str, T num);
  ///
  template<typename T> bool isEqual(cString &str, T num);
  
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

  // Wrapper for vector<vector<T>>
  template<typename T>
  class Vector2DContiguous {
  protected:
    const vector<T> v;
    const size_t s1;
    const size_t s2;
  public:
    Vector2DContiguous(size_t s1_, size_t s2_) : v(s1_*s2_,0), s1(s1_), s2(s2_) {;};
    Vector2DContiguous() : Vector2DContiguous(0,0) {;};
    size_t size(const size_t i) const { return (i == 0) ? s1 : s2; };
    vector<T>& operator()(const size_t i, const size_t j) { return v[i + j*s1]; };
    const vector<T>& operator()(const int i, const size_t j) const { return v[i + j*s1]; };
    // typename decltype(v)::iterator begin() { return v.begin(); };
    // typename decltype(v)::iterator end() { return v.end(); };
    // typename decltype(v)::const_iterator begin() const { return v.begin(); };
    // typename decltype(v)::const_iterator end() const { return v.end(); };
    // vector<T> flatten() const {
    //   vector<T> out;
    //   for (const auto &w : v) { out.insert(out.end(), w.begin(), w.end()); };
    //   return out;
    // }
  };
  
  // Wrapper for vector<vector<T>>
  template<typename T>
  class Vector2D {
  protected:
    vector<vector<T>> v;
  public:
    size_t size() const { return v.size(); };
    size_t size(const int i) const { return v[i].size(); };
    void resize(const int s1) { v.resize(s1); };
    void resize(const int s1, const int s2) {
      resize(s1);
      for_each(v.begin(), v.end(), [&](vector<T> &vT){ vT.resize(s2); });
    }
    vector<T>& operator[](const int i) { return v[i]; };
    const vector<T>& operator[](const int i) const { return v[i]; };
    typename decltype(v)::iterator begin() { return v.begin(); };
    typename decltype(v)::iterator end() { return v.end(); };
    typename decltype(v)::const_iterator begin() const { return v.begin(); };
    typename decltype(v)::const_iterator end() const { return v.end(); };
    vector<T> flatten() const {
      vector<T> out;
      for (const auto &w : v) { out.insert(out.end(), w.begin(), w.end()); };
      return out;
    }
  };
  
  // Wrapper for vector<vector<vector<T>>>
  template<typename T>
  class Vector3D {
  protected:
    vector<Vector2D<T>> v;
  public:
    size_t size() const { return v.size(); };
    size_t size(const int i) const { return v[i].size(); };
    size_t size(const int i, const int j) const { return v[i][j].size(); };
    void resize(const int s1) { v.resize(s1); };
    void resize(const int s1, const int s2) {
      resize(s1);
      for_each(v.begin(), v.end(), [&](vector<T> &vT){ vT.resize(s2); });
    }
    void resize(const int s1, const int s2, const int s3) {
      resize(s1);
      for_each(v.begin(), v.end(), [&](Vector2D<T> &vT){ vT.resize(s2, s3); });
    }
    Vector2D<T>& operator[](const int i) { return v[i]; };
    const Vector2D<T>& operator[](const int i) const { return v[i]; };
    typename decltype(v)::iterator begin() { return v.begin(); };
    typename decltype(v)::iterator end() { return v.end(); };
    typename decltype(v)::const_iterator begin() const { return v.begin(); };
    typename decltype(v)::const_iterator end() const { return v.end(); };
  };
    
}


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
    const Interpolator &ssfi;
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
		   const Interpolator &ssfi_,
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
    const Interpolator &ssfi;
    // Integrand
    double integrand(const double y) const ;
    // Compute static structure factor
    double ssf(double y_) const ;
  
  public:

    // Constructor
    Rdf(const double r_,
	const double cutoff_,
	const Interpolator &ssfi_,
	Integrator1DFourier &itg_) : r(r_), cutoff(cutoff_), itg(itg_),
				     ssfi(ssfi_)
    {

    };
    // Get result of integration 
    double get() const;
  
  };
  
}


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
