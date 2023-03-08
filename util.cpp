#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include "util.hpp"

namespace inpututil {

  vector<string> tokenize(cString &str, const char separator){
    stringstream strStream(str);
    string token;
    vector<string> tokens;
    while(getline(strStream, token, separator)) {
      tokens.push_back(token);
    }
    return tokens;
  }
  
  void matchKeyAndData(cVector<string> &keyword,
		       cString &input,
		       map<string,function<void(cString&, cString&)>> &funcArr){
    if (funcArr.find(keyword[0]) != funcArr.end()){
      funcArr[keyword[0]](keyword[1], input);
    }
    else {
      throw runtime_error("Unknown keyword: " + keyword[0] + "." + keyword[1]);
    }
  }
  
  void matchKeyAndData(cString &keyword,
		       cString &input,
		       map<string,function<void(cString&)>> &funcArr){
    if (funcArr.find(keyword) != funcArr.end()){
      funcArr[keyword](input);
    }
    else {
      throw runtime_error("Unknown keyword: " + keyword);
    }
  }

  template<> bool isNegative<int>(cString &str) {
    return stoi(str)<0;
  }

  template<> bool isNegative<double>(cString &str) {
    return stod(str)<0;
  }

  template<> bool isNotPositive<int>(cString &str) {
    return stoi(str)<=0;
  }

  template<> bool isNotPositive<double>(cString &str) {
    return stod(str)<=0;
  }

  template<> bool isLarger<int>(cString &str, int num){
    return stoi(str)>num;
  }

  template<> bool isLarger<double>(cString &str, double num){
    return stod(str)>num;
  }

  template<> bool isEqual(cString &str, int num){
    return stoi(str)==num;
  }
  
}


namespace vecUtil {

  // Element-wise sum between two vectors
  vector<double> sum(const vector<double> &v1,
		     const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res;
    transform(v1.begin(), v1.end(),
	      v2.begin(), back_inserter(res),
	      plus<double>());
    return res;
  }

  // Element-wise difference between two vectors
  vector<double> diff(const vector<double> &v1,
		      const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res;
    transform(v1.begin(), v1.end(),
	      v2.begin(), back_inserter(res),
	      minus<double>());
    return res;
  }
  
  // Root square difference between two vectors
  double rms(const vector<double> &v1,
	     const vector<double> &v2,
	     const bool normalize) {
    const vector<double> tmp = diff(v1,v2);
    double rms = inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.0);
    if (normalize) rms /= tmp.size();
    return sqrt(rms);
  }

  // Element-wise multiplication of a vector and a scalar
  vector<double> mult(const vector<double> &v,
		    const double a) {
    vector<double> res = v;
    transform(res.begin(), res.end(), res.begin(), [&a](double c){return c*a;});
    return res;
  }

  // Vector2D class
  size_t Vector2D::size() const {
    return s1*s2;
  }
  
  size_t Vector2D::size(const size_t i) const {
    assert(i == 0 || i == 1);
    return (i == 0) ? s1 : s2;
  }

  void Vector2D::resize(const size_t s1_, const size_t s2_) {
    if (s1_ == s1 && s2_ == s2) { return; }
    v.clear();
    v.resize(s1_*s2_, 0.0);
  }

  double& Vector2D::operator()(const size_t i, const size_t j) {
    return v[j + i*s2];
  }
  
  const double& Vector2D::operator()(const size_t i, const size_t j) const {
    return v[j + i*s2];
  }

  const double& Vector2D::operator()(const size_t i) const {
    return v[i*s2];
  }

  vector<double>::iterator Vector2D::begin() {
    return v.begin();
  }

  vector<double>::iterator Vector2D::end() {
    return v.end();
  }

  vector<double>::const_iterator Vector2D::begin() const {
    return v.begin();
  }

  vector<double>::const_iterator Vector2D::end() const {
    return v.end();
  }
  
  
  void Vector2D::fill(const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi){ vi = num;});
  }

  void Vector2D::fillRow(const size_t i, const double &num) {
    const auto &dest = v.begin() + i*s2;
    std::for_each(dest, dest + s2, [&](double &vi){ vi = num;});
  }
  
  void Vector2D::fillRow(const size_t i, const vector<double> &num) {
    assert(num.size() == s2);
    std::copy(num.begin(), num.end(), v.begin() + i*s2);
  }
  
  
}

namespace thermoUtil {

  double InternalEnergy::ssf(double y) const {
    return ssfi.eval(y);
  }
  
  double InternalEnergy::integrand(double y) const {
    return ssf(y) - 1;
  }

  double InternalEnergy::get() const  {
    auto func = [&](double y)->double{return integrand(y);};
    itg.compute(func, yMin, yMax);
    return itg.getSolution()/(M_PI * rs * lambda);
  }

  double Rdf::ssf(double y) const {
    return ssfi.eval(y);
  }

  double Rdf::integrand(double y) const {
    if (y > cutoff) return 0;
    return y*(ssf(y) - 1);
  }

  double Rdf::get() const {
    assert(r>0);
    auto func = [&](double y)->double{return integrand(y);};
    itg.setR(r);
    itg.compute(func);
    return 1 + 1.5 * itg.getSolution()/r;
  }

}
