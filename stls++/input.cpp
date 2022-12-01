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
  
}

void Input::assignInputToStlsData(const string &keyword, const string &value){
  cout << keyword << " = " << value << endl;
}

// --- staticInput ---

staticInput::staticInput(){
  aMix = 1.0;
  errMin = 1e-5;
  dx = 0.1;
  xmax = 10.0;
  vector<double> muGuessDefault = {-10, 10};
  muGuess.assign(muGuessDefault.begin(), muGuessDefault.end());
  nl = 128;
  nIter = 1000;
}

double staticInput::getMixingParameter(){
  return aMix;
}

double staticInput::getErrMin(){
  return errMin;
}

double staticInput::getWaveVectorGridRes(){
  return dx;
}

double staticInput::getWaveVectorGridCutoff(){
  return xmax;
}
 
vector<double> staticInput::getChemicalPotentialGuess(){
  return muGuess;
}

size_t staticInput::getNMatsubara(){
  return nl;
}

size_t staticInput::getNIter(){
  return nIter;
}

void staticInput::setMixingParameter(double aMix){
  this->aMix = aMix;
}

void staticInput::setErrMin(double errMin){
  this->errMin = errMin;
}

void staticInput::setChemicalPotentialGuess(vector<double> muGuess){
  this->muGuess.assign(muGuess.begin(), muGuess.end());
}
 
void staticInput::setWaveVectorGridRes(double dx){
  this->dx = dx;
}

void staticInput::setWaveVectorGridCutoff(double xmax){
  this->xmax = xmax;
}

void staticInput::setNMatsubara(size_t nl){
  this->nl = nl;
}

void staticInput::setNIter(size_t nIter){
  this->nIter = nIter; 
}

// --- stlsInput ---

stlsInput::stlsInput(){
  IETMapping = "standard";
  restartFileName = NO_FILE_NAME;
}

string stlsInput::getIETMapping(){
  return IETMapping;
}

string stlsInput::getRestartFileName(){
  return restartFileName;
}

void stlsInput::setIETMapping(string IETMapping){
  this->IETMapping = IETMapping;
}

void stlsInput::setRestartFileName(string restartFileName){
  this->restartFileName = restartFileName;
} 
  
 
