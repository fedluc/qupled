#include <omp.h>
#include "esa.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

// -----------------------------------------------------------------
// ESA class
// -----------------------------------------------------------------

// ESA compute function
int ESA::compute(){
  try {
    init();  
    if (verbose) cout << "Structural properties calculation ..." << endl;
    doESA();
    if (verbose) cout << "Done" << endl;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}


// ESA main function
void ESA::doESA() {

// Start timing
double tic = omp_get_wtime();
// Compute ESA static local field correction
computeSlfc();
// Compute static structure factor
computeSsf();
// End timing
double toc = omp_get_wtime();
// Print diagnostic
if (verbose) {
    printf("Elapsed time: %f seconds\n", toc - tic);
    fflush(stdout);
}
  
}


void ESA::computeSlfc() {

    // Get the degeneracy and coupling parameters
    const double theta = in.getDegeneracy();
    const double rs = in.getCoupling();

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

    const double theta_pow = pow(theta, 3.0/2.0);
    const double aa = aa1 + aa2 * theta + aa3 * theta_pow;
    const double ba = ba1 + ba2 * theta + ba3 * theta_pow;
    const double ca = ca1 + ca2 * theta + ca3 * theta_pow;
    const double ab = ab1 + ab2 * theta + ab3 * theta_pow;
    const double bb = bb1 + bb2 * theta + bb3 * theta_pow;
    const double cb = cb1 + cb2 * theta + cb3 * theta_pow;
    const double ac = ac1 + ac2 * theta + ac3 * theta_pow;
    const double bc = bc1 + bc2 * theta + bc3 * theta_pow;
    const double cc = cc1 + cc2 * theta + cc3 * theta_pow;
    const double ad = ad1 + ad2 * theta + ad3 * theta_pow;
    const double bd = bd1 + bd2 * theta + bd3 * theta_pow;
    const double cd = cd1 + cd2 * theta + cd3 * theta_pow;
    const double xm = Ax + Bx * theta + Cx * pow(theta, 2.0); 

    const double g0a = (g0a0 + g0aa1 * theta)/(1.0 + g0ba1 * theta + g0ba2 * pow(theta, 3.0));
    const double g0b = (g0b0 + g0ab1 * sqrt(theta))/(1.0 + g0bb1 * theta + g0bb2 * pow(theta, 2.0));
    const double g0c = (g0c0 + g0ac1 * sqrt(theta) + g0ac2 * theta_pow)/(1.0 + g0bc1 * theta + g0bc2 * pow(theta, 2.0));
    const double g0d = (g0d0 + g0ad1 * sqrt(theta))/(1.0 + g0bd1 * theta + g0bd2 * pow(theta, 2.0));
    const double g = 1.0/2.0 * (1.0 + g0a * sqrt(rs) + g0b * rs)/(1.0 + g0c * rs + g0d * pow(rs, 3.0));
    
    const double a = (aa + ba * rs)/(1.0 + ca * rs);
    const double b = (ab + bb * rs)/(1.0 + cb * rs);
    const double c = (ac + bc * rs)/(1.0 + cc * rs);
    const double d = (ad + bd * rs)/(1.0 + cd * rs);

    // Assigning derivatives to variables
    const double dfxc_rs = (fxc(theta, rs + dx) - fxc(theta, rs - dx)) / (2.0 * dx);
    const double dfxc_t = (fxc(theta + dx, rs) - fxc(theta - dx, rs)) / (2.0 * dx);
    const double dfxc_rs2 = (fxc(theta, rs + dx) - 2.0 * fxc(theta, rs) + fxc(theta, rs - dx)) / (dx * dx);
    const double dfxc_t2 = (fxc(theta + dx, rs) - 2.0 * fxc(theta, rs) + fxc(theta - dx, rs)) / (dx * dx);
    const double dfxc_t_rs = (fxc(theta + dx, rs + dx) - fxc(theta + dx, rs - dx) - fxc(theta - dx, rs + dx) + fxc(theta - dx, rs - dx)) / (4.0 * dx * dx);
    
    const int nx = wvg.size();
    slfc.resize(nx);
    // Loop over the wave vector grid size
    for (int i=0; i<nx; ++i) {
                
        const double AF = 1.0/2.0 * (1.0 + tanh(Eta * (wvg[i] - xm)));

        const double GCSR = -(M_PI/12.0) * lambda * rs * pow(wvg[i], 2.0) * (4.0 * pow(theta, 2.0) * dfxc_t2 + pow(rs, 2.0) * dfxc_rs2 + 4.0 * theta * rs * dfxc_t_rs - 2.0 * theta * dfxc_t - 2.0 * rs * dfxc_rs);
                
        const double Gnnfit = (1.0 + a * wvg[i] + b * pow(wvg[i], 1.0/2.0))/(1.0 + c * wvg[i] + d * pow(wvg[i], 1.25) + GCSR);
        
        const double Gesa = GCSR * Gnnfit * (1.0 - AF) + (1.0 - g) * AF;
        
        slfc[i] = Gesa;
    }
    // Assign slfc values to the slfcOld vector
    slfcOld = slfc;
}

// QMC free energy function constants
double ESA::fxc(double theta, double rs) const {

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

    double fa = 0.610887 * tanh(pow(theta, -1.0)) * (0.75 + 3.04363 * pow(theta, 2.0) - 0.09227 * pow(theta, 3.0) + 1.7035 * pow(theta, 4.0)) / (1.0 + 8.31051 * pow(theta, 2.0) + 5.1105 * pow(theta, 4.0));
    double fb = tanh(1.0 / sqrt(theta)) * (fb1 + fb2 * pow(theta, 2.0) + fb3 * pow(theta, 4.0)) / (1.0 + fb4 * pow(theta, 2.0) + fb5 * pow(theta, 4.0));
    double fd = tanh(1.0 / sqrt(theta)) * (fd1 + fd2 * pow(theta, 2.0) + fd3 * pow(theta, 4.0)) / (1.0 + fd4 * pow(theta, 2.0) + fd5 * pow(theta, 4.0));
    double fe = tanh(1.0 / theta) * (fe1 + fe2 * pow(theta, 2.0) + fe3 * pow(theta, 4.0)) / (1.0 + fe4 * pow(theta, 2.0) + fe5 * pow(theta, 4.0));
    double fc = (fc1 + fc2 * exp(-(1.0 / theta))) * fe;

    return -(1.0 / rs) * (omega * fa + fb * sqrt(rs) + fc * rs) / (1.0 + fd * sqrt(rs) + fe * rs);
    }
