#include <string>
#include <fstream>
#include <sstream>
#include "input.hpp"

// --- Input ---
 
Input::Input(){
  theory = "stls";
  Theta = 1.0;
  rs = 1.0;
  nThreads = 1;
}

string Input::getTheory(){
  return theory;
}

double Input::getDegeneracy(){
  return Theta;
}

double Input::getCoupling(){
  return rs;
}

void Input::setTheory(cString &theory){
  if (theory != "stls") {
    throw runtime_error("Unknown theory: " + theory);
  }
  this->theory = theory;
}

void Input::setDegeneracy(cString &Theta){
  double ThetaNum = stod(Theta);
  if (ThetaNum < 0.0) {
    throw runtime_error("The quantum degeneracy parameter can't be negative");
  }
  this->Theta = ThetaNum;
}

void Input::setCoupling(cString &rs){
  double rsNum = stod(rs);
  if (rsNum <= 0.0) {
    throw runtime_error("The quantum coupling parameter must be larger than zero");
  }
  this->rs = rsNum;
}

void Input::setThreads(cString  &nThreads){
  int nThreadsNum = stoi(nThreads);
  if (nThreadsNum <= 0.0) {
    throw runtime_error("The number of threads must be positive");
  }
  this->nThreads = nThreadsNum;
}

void Input::readInput(cString &fileName){
  ifstream file(fileName);
  if (file.is_open()) {
    string line;
    while (getline(file, line)) {
      parseInputLine(line);
    }
    file.close();
  }
  else {
    throw runtime_error("Input file " + fileName + " could not be opened.");    
  }   
}

void Input::parseInputLine(cString &line){
  bool isComment = line[0] == '#';
  bool isEmpty = line[0] == '\n' || line.length()==0;
  if (!isComment && !isEmpty) {
    vector<string> tokens = tokenize(line, ' ');
    if (tokens.size() < 2) {
      throw runtime_error("wrong line format: " + line);
    }
    assignInputToData(tokens);    
  }
}

vector<string> Input::tokenize(cString &str, const char separator){
 stringstream strStream(str);
 string token;
 vector<string> tokens;
 while(getline(strStream, token, separator)) {
   tokens.push_back(token);
 }
 return tokens;
}

void Input::assignInputToData(cVector<string> &input){
  cVector<string> keyword = tokenize(input[0], '.');
  map<string, function<void(cString &, cString &)>> funcArr;
  funcArr["base"] = [this](cString &s1, cString &s2) {this->assignInputToBaseData(s1, s2);};
  funcArr["static"] = [this](cString &s1, cString &s2) {this->assignInputToStaticData(s1, s2);};
  funcArr["stls"] = [this](cString &s1, cString &s2) {this->assignInputToStlsData(s1, s2);};
  try{
    matchKeyAndData(keyword, input, funcArr);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
}

void Input::assignInputToBaseData(cString &keyword, cString &value){
  map<string, function<void(cString &)>> funcArr;
  funcArr["theory"] = [this](cString &s1) {this->setTheory(s1);};
  funcArr["degeneracy"] = [this](cString &s1) {this->setDegeneracy(s1);};
  funcArr["coupling"] = [this](cString &s1) {this->setCoupling(s1);};
  funcArr["threads"] = [this](cString &s1) {this->setThreads(s1);};
  try{
    matchKeyAndData(keyword, value, funcArr);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
}

void Input::assignInputToStaticData(cString &keyword, cString &value){
  stat.assignInputToData(keyword, value);
}

void Input::assignInputToStlsData(cString &keyword, cString &value){
  cout << keyword << " = " << value << endl;
}

void Input::matchKeyAndData(cVector<string> &keyword,
			    cVector<string> &input,
			    map<string,function<void(cString&, cString&)>> &funcArr){
  if (funcArr.find(keyword[0]) != funcArr.end()){
    funcArr[keyword[0]](keyword[1], input[1]);
  }
  else {
    throw runtime_error("Unknown keyword: " + input[0]);
  }
}

void Input::matchKeyAndData(cString &keyword,
			    cString &input,
			    map<string,function<void(cString&)>> &funcArr){
  if (funcArr.find(keyword) != funcArr.end()){
    funcArr[keyword](input);
  }
  else {
    throw runtime_error("Unknown keyword: " + keyword);
  }
}

// --- StaticInput ---

StaticInput::StaticInput(){
  aMix = 1.0;
  errMin = 1e-5;
  dx = 0.1;
  xmax = 10.0;
  vector<double> muGuessDefault = {-10, 10};
  muGuess.assign(muGuessDefault.begin(), muGuessDefault.end());
  nl = 128;
  nIter = 1000;
}

double StaticInput::getMixingParameter(){
  return aMix;
}

double StaticInput::getErrMin(){
  return errMin;
}

double StaticInput::getWaveVectorGridRes(){
  return dx;
}

double StaticInput::getWaveVectorGridCutoff(){
  return xmax;
}
 
vector<double> StaticInput::getChemicalPotentialGuess(){
  return muGuess;
}

size_t StaticInput::getNMatsubara(){
  return nl;
}

size_t StaticInput::getNIter(){
  return nIter;
}

void StaticInput::setMixingParameter(const double aMix){
  if (aMix < 0.0 || aMix > 1.0) {
    throw runtime_error("The mixing parameter must be a number between zero and one");
  }
  this->aMix = aMix;
}

void StaticInput::setErrMin(const double errMin){
  if (errMin <= 0.0) {
    throw runtime_error("The minimum error for convergence must be larger than zero");
  }
  this->errMin = errMin;
}

void StaticInput::setChemicalPotentialGuess(const string &muGuessStr){
  vector<string> muGuessVecStr = Input::tokenize(muGuessStr, ',');
  if (muGuessVecStr.size() != 2) {
    throw runtime_error("Wrong format for the chemical potential input.");
  }
  vector<double> muGuess;
  for (string mu : muGuessVecStr) muGuess.push_back(stod(mu));
  this->muGuess.assign(muGuess.begin(), muGuess.end());
}
 
void StaticInput::setWaveVectorGridRes(const double dx){
  if (dx <= 0.0) {
    throw runtime_error("The wave-vector grid resolution must be larger than zero");
  }
  this->dx = dx;
}

void StaticInput::setWaveVectorGridCutoff(const double xmax){
  if (xmax <= 0.0) {
    throw runtime_error("The wave-vector grid cutoff must be larger than zero");
  }
  if (xmax < dx) {
    throw runtime_error("The wave-vector grid cutoff must be larger than the resolution");
  }
  this->xmax = xmax;
}

void StaticInput::setNMatsubara(const size_t nl){
  if (nl < 0.0) {
    throw runtime_error("The number of matsubara frequencies can't be negative");
  }
  this->nl = nl;
}

void StaticInput::setNIter(const size_t nIter){
  if (nIter < 0.0) {
    throw runtime_error("The maximum number of iterations can't be negative");
  }
  this->nIter = nIter; 
}

void StaticInput::assignInputToData(const string &keyword, const string &value){
  // if (keyword == allowedKeywords[0])
  //   setMixingParameter(stod(value));
  // else if (keyword == allowedKeywords[1])
  //   setErrMin(stod(value));
  // else if (keyword == allowedKeywords[2])
  //   setWaveVectorGridRes(stod(value));
  // else if (keyword == allowedKeywords[3])
  //   setWaveVectorGridCutoff(stod(value));
  // else if (keyword == allowedKeywords[4])
  //   setChemicalPotentialGuess(value);
  // else if (keyword == allowedKeywords[5])
  //   setNMatsubara(stoi(value));
  // else if (keyword == allowedKeywords[6])
  //   setNIter(stoi(value));
  // else
  //   throw runtime_error("Unknown keyword: " + keyword);
}

// --- StlsInput ---

StlsInput::StlsInput(){
  IETMapping = "standard";
  restartFileName = NO_FILE_NAME;
}

string StlsInput::getIETMapping(){
  return IETMapping;
}

string StlsInput::getRestartFileName(){
  return restartFileName;
}

void StlsInput::setIETMapping(string IETMapping){
  this->IETMapping = IETMapping;
}

void StlsInput::setRestartFileName(string restartFileName){
  this->restartFileName = restartFileName;
} 
  
 
