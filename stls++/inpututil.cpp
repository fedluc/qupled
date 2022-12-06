#include <string>
#include <fstream>
#include <sstream>
#include "inpututil.hpp"

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

  template <> bool isNegative<int>(cString &str) {
    return stoi(str)<0;
  }

  template <> bool isNegative<double>(cString &str) {
    return stod(str)<0;
  }

  template <> bool isNotPositive<int>(cString &str) {
    return stoi(str)<=0;
  }

  template <> bool isNotPositive<double>(cString &str) {
    return stod(str)<=0;
  }

  template <> bool isLarger<int>(cString &str1, int num){
    return stoi(str1)>num;
  }

  template <> bool isLarger<double>(cString &str1, double num){
    return stod(str1)>num;
  }
  
}
