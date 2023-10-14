#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------


class StructProp {
  // Legend : rsm = rs - drs, rsp = rs + drs,
  // tm = theta - dtheta, tp = theta + dtheta 
private:
  // Stls rsmtm;
  // Stls rstm;
  // Stls rsptmp;
  // Stls rsmt;
  // Stls rst;
  // Stls rspt;
  // Stls rsmtp;
  // Stls rstp;
  // Stls rsptp;
public:
  StructProp(const StlsInput &in_);
};

class ThermoProp {
private:
  std::vector<double> rstm;
  std::vector<double> rts;
  std::vector<double> rstp;
};

class VSStls {

private: 


  // Private members
  StructProp stls;
  ThermoProp rsu;
  ThermoProp rsa;
  
    
public:

  // Constructors
  VSStls(const StlsInput &in_) : stls(in), rsu(), rsa() { ; };
  // Compute stls scheme
  int compute();
  
};


#endif
