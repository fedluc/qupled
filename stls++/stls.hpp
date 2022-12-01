#include <vector>

#ifndef STLS_HPP
#define STLS_HPP

using namespace std;

class stls {

private: 

  // Wave vector grid
  vector<double> wvg;
  // Ideal density response
  vector<double> idr;
  // Static local field correction
  vector<double> slfc;
  // Static structure factor
  vector<double> ssf;
  // Hartree-Fock static structure factor
  vector<double> ssfHF;
 
public:

  // Constructor
  stls();



};

#endif