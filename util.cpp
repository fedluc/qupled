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
    itg->compute(func, yMin, yMax);
    return itg->getSolution()/(M_PI * rs * lambda);
  }

  double Rdf::ssf(double y) const {
    return ssfi.eval(y);
  }

  double Rdf::integrand(double y) const {
    if (y > cutoff) return 0;
    return y*(ssf(y) - 1);
  }

  double Rdf::get() const {
    auto func = [&](double y)->double{return integrand(y);};
    itg->compute(func);
    return 1 + 1.5 * itg->getSolution()/r;
  }

}
