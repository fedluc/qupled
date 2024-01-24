#include <omp.h>
#include "esa.hpp"

using namespace std;

// -----------------------------------------------------------------
// ESA class
// -----------------------------------------------------------------

ESA::ESA(const RpaInput &in_) : Rpa(in_, true, false) {
  if (verbose) cout << "Computing static local field correction: ";
  slfc.resize(wvg.size());
  computeSlfc();
  if (verbose) cout << "Done" << endl;
  if (verbose) cout << "Computing static structure factor: "; 
  computeSsf();
  if (verbose) cout << "Done" << endl;
}

void ESA::computeSlfc() {

    // Get the degeneracy and coupling parameters
    const double theta = in.getDegeneracy();
    const double rs = in.getCoupling();
    const double theta3_2 = pow(theta, 1.5);
    const double theta2 = theta * theta;
    const double theta3 = theta * theta2;
    const double sqrtTheta = sqrt(theta);
    const double sqrtRs = sqrt(rs);
    const double rs2 = rs * rs;
    const double rs3 = rs2 * rs;
    const double lambdaRs = -(M_PI/12.0) * lambda * rs;
    
    // Constants for the computation of the static local field correction
    constexpr double g0aa1 = 18.4376509088802;
    constexpr double g0ba1 = 24.1338558554951;
    constexpr double g0ba2 = 1.86499223116244;
    constexpr double g0a0 = 0.18315;
    constexpr double g0ab1 = -0.243679875065179;
    constexpr double g0bb1 = 0.252577;
    constexpr double g0bb2 = 0.12704315679703;
    constexpr double g0b0 = -0.0784043;
    constexpr double g0ac1 = 2.23662699250646;
    constexpr double g0ac2 = 0.448936594134834;
    constexpr double g0bc1 = 0.445525554738623;
    constexpr double g0bc2 = 0.408503601408896;
    constexpr double g0c0 = 1.02232;
    constexpr double g0ad1 = 0.0589015402409623;
    constexpr double g0bd1 = -0.598507865444083;
    constexpr double g0bd2 = 0.513161640681146;
    constexpr double g0d0 = 0.0837741;

    constexpr double Eta = 3.0;
    constexpr double Ax = 2.64;
    constexpr double Bx = 0.31;
    constexpr double Cx = 0.08;
    constexpr double aa1 = 0.66477593;
    constexpr double aa2 = -4.59280227;
    constexpr double aa3 = 1.24649624;
    constexpr double ba1 = -1.27089927;
    constexpr double ba2 = 1.26706839;
    constexpr double ba3 = -0.4327608;
    constexpr double ca1 = 2.09717766;
    constexpr double ca2 = 1.15424724;
    constexpr double ca3 = -0.65356955;
    constexpr double ab1 = -1.0206202;
    constexpr double ab2 = 5.16041218;
    constexpr double ab3 = -0.23880981;
    constexpr double bb1 = 1.07356921;
    constexpr double bb2 = -1.67311761;
    constexpr double bb3 = 0.58928105;
    constexpr double cb1 = 0.8469662;
    constexpr double cb2 = 1.54029035;
    constexpr double cb3 = -0.71145445;
    constexpr double ac1 = -2.31252076;
    constexpr double ac2 = 5.83181391;
    constexpr double ac3 = 2.29489749;
    constexpr double bc1 = 1.76614589;
    constexpr double bc2 = -0.09710839;
    constexpr double bc3 = -0.33180686;
    constexpr double cc1 = 0.56560236;
    constexpr double cc2 = 1.10948188;
    constexpr double cc3 = -0.43213648;
    constexpr double ad1 = 1.3742155;
    constexpr double ad2 = -4.01393906;
    constexpr double ad3 = -1.65187145;
    constexpr double bd1 = -1.75381153;
    constexpr double bd2 = -1.17022854;
    constexpr double bd3 = 0.76772906;
    constexpr double cd1 = 0.63867766;
    constexpr double cd2 = 1.07863273;
    constexpr double cd3 = -0.35630091;

    const double aa = aa1 + aa2 * theta + aa3 * theta3_2;
    const double ba = ba1 + ba2 * theta + ba3 * theta3_2;
    const double ca = ca1 + ca2 * theta + ca3 * theta3_2;
    const double ab = ab1 + ab2 * theta + ab3 * theta3_2;
    const double bb = bb1 + bb2 * theta + bb3 * theta3_2;
    const double cb = cb1 + cb2 * theta + cb3 * theta3_2;
    const double ac = ac1 + ac2 * theta + ac3 * theta3_2;
    const double bc = bc1 + bc2 * theta + bc3 * theta3_2;
    const double cc = cc1 + cc2 * theta + cc3 * theta3_2;
    const double ad = ad1 + ad2 * theta + ad3 * theta3_2;
    const double bd = bd1 + bd2 * theta + bd3 * theta3_2;
    const double cd = cd1 + cd2 * theta + cd3 * theta3_2;
    const double xm = Ax + Bx * theta + Cx * theta2; 

    const double g0a = (g0a0 + g0aa1 * theta)/(1.0 + g0ba1 * theta + g0ba2 * theta3);
    const double g0b = (g0b0 + g0ab1 * sqrtTheta)/(1.0 + g0bb1 * theta + g0bb2 * theta2);
    const double g0c = ( (g0c0 + g0ac1 * sqrtTheta + g0ac2 * theta3_2)
			 /(1.0 + g0bc1 * theta + g0bc2 * theta2) );
    const double g0d = (g0d0 + g0ad1 * sqrtTheta)/(1.0 + g0bd1 * theta + g0bd2 * theta2);
    const double g = 1.0/2.0 * (1.0 + g0a * sqrtRs + g0b * rs)/(1.0 + g0c * rs + g0d * rs3);
    
    const double a = (aa + ba * rs)/(1.0 + ca * rs);
    const double b = (ab + bb * rs)/(1.0 + cb * rs);
    const double c = (ac + bc * rs)/(1.0 + cc * rs);
    const double d = (ad + bd * rs)/(1.0 + cd * rs);

    // Assigning derivatives to variables
    vector<double> fxcData;
    for (const double& thetaTmp : {theta - dx, theta, theta + dx}) {
      for (const double& rsTmp : {rs - dx, rs, rs + dx}){
	fxcData.push_back(fxc(thetaTmp, rsTmp));
      }
    }
    const double tdx = 2.0 * dx;
    const double dx2 = dx * dx;
    const double fdx2 = 4.0 * dx2;
    const double fxcr = (fxcData[5] - fxcData[3]) / tdx;
    const double fxct = (fxcData[7] - fxcData[1]) / tdx;
    const double fxcrr = (fxcData[5] - 2.0 * fxcData[4] + fxcData[3]) / dx2;
    const double fxctt = (fxcData[7] - 2.0 * fxcData[4] + fxcData[1]) / dx2;
    const double fxctr = (fxcData[8] - fxcData[6] - fxcData[2] + fxcData[0]) / fdx2;
    
    // Loop over the wave vector grid size
    for (size_t i = 0; i < wvg.size(); ++i) {
      const double& x = wvg[i];
      const double x2 = x * x;
      const double sqrtX = sqrt(x);
      const double xp125 = pow(x, 1.25);
      const double AF = 1.0/2.0 * (1.0 + tanh(Eta * (x - xm)));
      const double GCSR = lambdaRs * x2 * (4.0 * theta2 * fxctt + rs2 * fxcrr
					   + 4.0 * theta * rs * fxctr
					   - 2.0 * theta * fxct - 2.0 * rs * fxcr);
      const double Gnnfit = (1.0 + a * x + b * sqrtX)/(1.0 + c * x + d * xp125 + GCSR);
      const double Gesa = GCSR * Gnnfit * (1.0 - AF) + (1.0 - g) * AF;
      slfc[i] = Gesa;
    }
}

// QMC free energy function constants
double ESA::fxc(const double& theta,
		const double& rs) const {

  const double thetaInv = 1.0/theta;
  const double theta2 = theta * theta;
  const double theta3 = theta * theta2;
  const double theta4 = theta2 * theta2;
  const double tanhThetaInv = tanh(thetaInv);
  const double tanhSqrtThetaInv = tanh(sqrt(thetaInv));
  const double rsInv = 1.0/rs;
  const double sqrtRs = sqrt(rs);
  constexpr double omega = 1.0;
  constexpr double fb1 = 0.3436902;
  constexpr double fb2 = 7.82159531356;
  constexpr double fb3 = 0.300483986662;
  constexpr double fb4 = 15.8443467125;
  const double fb5 = fb3*pow(3.0/2.0, 1.0/2.0)*omega/lambda;
  constexpr double fc1 = 0.8759442;
  constexpr double fc2 = -0.230130843551;
  constexpr double fd1 = 0.72700876;
  constexpr double fd2 = 2.38264734144;
  constexpr double fd3 = 0.30221237251;
  constexpr double fd4 = 4.39347718395;
  constexpr double fd5 = 0.729951339845;
  constexpr double fe1 = 0.25388214;
  constexpr double fe2 = 0.815795138599;
  constexpr double fe3 = 0.0646844410481;
  constexpr double fe4 = 15.0984620477;
  constexpr double fe5 = 0.230761357474;

  const double fa = 0.610887 * tanhThetaInv * ( (0.75 + 3.04363 * theta2 - 0.09227 * theta3 + 1.7035 * theta4)
						/ (1.0 + 8.31051 * theta2 + 5.1105 * theta4) );
  const double fb = tanhSqrtThetaInv * ( (fb1 + fb2 * theta2 + fb3 * theta4)
					 / (1.0 + fb4 * theta2 + fb5 * theta4) );
  const double fd = tanhSqrtThetaInv * ( (fd1 + fd2 * theta2 + fd3 * theta4)
					 / (1.0 + fd4 * theta2 + fd5 * theta4) );
  const double fe = tanhThetaInv * ( (fe1 + fe2 * theta2 + fe3 * theta4)
				     / (1.0 + fe4 * theta2 + fe5 * theta4) );
  const double fc = (fc1 + fc2 * exp(-thetaInv)) * fe;

  return -rsInv * (omega * fa + fb * sqrtRs + fc * rs) / (1.0 + fd * sqrtRs + fe * rs);
    
}
