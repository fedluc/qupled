#ifndef INPUTUTIL_HPP
#define INPUTUTIL_HPP

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <functional>

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

#endif
