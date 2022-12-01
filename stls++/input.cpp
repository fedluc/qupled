#include <string>
#include <fstream>
#include <sstream>
#include <vector>
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

void Input::setTheory(const string &theory){
  if (theory != "stls") {
    throw runtime_error("Unknown theory: " + theory);
  }
  this->theory = theory;
}

void Input::setDegeneracy(const double Theta){
  if (Theta < 0.0) {
    throw runtime_error("The quantum degeneracy parameter can't be negative");
  }
  this->Theta = Theta;
}

void Input::setCoupling(const double rs){
  if (rs <= 0.0) {
    throw runtime_error("The quantum coupling parameter must be larger than zero");
  }
  this->rs = rs;
}

void Input::setThreads(const int nThreads){
  if (nThreads <= 0.0) {
    throw runtime_error("The number of threads must be positive");
  }
  this->nThreads = nThreads;
}

void Input::readInput(const string &fileName){
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

void Input::parseInputLine(const string &line){
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

vector<string> Input::tokenize(const string &str, const char separator){
 stringstream strStream(str);
 string token;
 vector<string> tokens;
 while(getline(strStream, token, separator)) {
   tokens.push_back(token);
 }
 return tokens;
}

void Input::assignInputToData(const vector<string> &input){
  const vector<string> keyword = tokenize(input[0], '.');
  try{
    if (keyword[0] == allowedKeywords[0])
      assignInputToBaseData(keyword[1], input[1]);
    else if (keyword[0] == allowedKeywords[5])
      assignInputToStaticData(keyword[1], input[1]);
    else if (keyword[0] == allowedKeywords[6])
      assignInputToStlsData(keyword[1], input[1]);
    else
      throw runtime_error("Unknown keyword: " + input[0]);  
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
  }
}

void Input::assignInputToBaseData(const string &keyword, const string &value){
  if (keyword == allowedKeywords[1])
    setTheory(value);
  else if (keyword == allowedKeywords[2])
    setDegeneracy(stod(value));
  else if (keyword == allowedKeywords[3])
    setCoupling(stod(value));
  else if (keyword == allowedKeywords[4])
    setThreads(stoi(value));
  else
    throw runtime_error("Unknown keyword: " + keyword);
}

void Input::assignInputToStaticData(const string &keyword, const string &value){
  stat.assignInputToData(keyword, value);
}

void Input::assignInputToStlsData(const string &keyword, const string &value){
  cout << keyword << " = " << value << endl;
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
  if (keyword == allowedKeywords[0])
    setMixingParameter(stod(value));
  else if (keyword == allowedKeywords[1])
    setErrMin(stod(value));
  else if (keyword == allowedKeywords[2])
    setWaveVectorGridRes(stod(value));
  else if (keyword == allowedKeywords[3])
    setWaveVectorGridCutoff(stod(value));
  else if (keyword == allowedKeywords[4])
    setChemicalPotentialGuess(value);
  else if (keyword == allowedKeywords[5])
    setNMatsubara(stoi(value));
  else if (keyword == allowedKeywords[6])
    setNIter(stoi(value));
  else
    throw runtime_error("Unknown keyword: " + keyword);
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
  
 
