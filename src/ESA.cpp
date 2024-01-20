#include <omp.h>
#include "ESA.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

// -----------------------------------------------------------------
// ESA class
// -----------------------------------------------------------------

ESA::ESA(const StlsInput& input) : Stls(input) {
   
}

// ESA compute function
int ESA::compute(){
  try {
    init();  
    if (verbose) cout << "Structural properties calculation ESA ..." << endl;
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
    const double g0aa1 = 18.4376509088802;
    const double g0ba1 = 24.1338558554951;
    const double g0ba2 = 1.86499223116244;
    const double g0a0 = 0.18315;
    const double g0ab1 = -0.243679875065179;
    const double g0bb1 = 0.252577;
    const double g0bb2 = 0.12704315679703;
    const double g0b0 = -0.0784043;
    const double g0ac1 = 2.23662699250646;
    const double g0ac2 = 0.448936594134834;
    const double g0bc1 = 0.445525554738623;
    const double g0bc2 = 0.408503601408896;
    const double g0c0 = 1.02232;
    const double g0ad1 = 0.0589015402409623;
    const double g0bd1 = -0.598507865444083;
    const double g0bd2 = 0.513161640681146;
    const double g0d0 = 0.0837741;

    const double Eta = 3;
    const double Ax = 2.64;
    const double Bx = 0.31;
    const double Cx = 0.08;
    const double aa1 = 0.66477593;
    const double aa2 = -4.59280227;
    const double aa3 = 1.24649624;
    const double ba1 = -1.27089927;
    const double ba2 = 1.26706839;
    const double ba3 = -0.4327608;
    const double ca1 = 2.09717766;
    const double ca2 = 1.15424724;
    const double ca3 = -0.65356955;
    const double ab1 = -1.0206202;
    const double ab2 = 5.16041218;
    const double ab3 = -0.23880981;
    const double bb1 = 1.07356921;
    const double bb2 = -1.67311761;
    const double bb3 = 0.58928105;
    const double cb1 = 0.8469662;
    const double cb2 = 1.54029035;
    const double cb3 = -0.71145445;
    const double ac1 = -2.31252076;
    const double ac2 = 5.83181391;
    const double ac3 = 2.29489749;
    const double bc1 = 1.76614589;
    const double bc2 = -0.09710839;
    const double bc3 = -0.33180686;
    const double cc1 = 0.56560236;
    const double cc2 = 1.10948188;
    const double cc3 = -0.43213648;
    const double ad1 = 1.3742155;
    const double ad2 = -4.01393906;
    const double ad3 = -1.65187145;
    const double bd1 = -1.75381153;
    const double bd2 = -1.17022854;
    const double bd3 = 0.76772906;
    const double cd1 = 0.63867766;
    const double cd2 = 1.07863273;
    const double cd3 = -0.35630091;

    double aa = aa1 + aa2 * theta + aa3 * pow(theta, 3.0/2.0);
    double ba = ba1 + ba2 * theta + ba3 * pow(theta, 3.0/2.0);
    double ca = ca1 + ca2 * theta + ca3 * pow(theta, 3.0/2.0);
    double ab = ab1 + ab2 * theta + ab3 * pow(theta, 3.0/2.0);
    double bb = bb1 + bb2 * theta + bb3 * pow(theta, 3.0/2.0);
    double cb = cb1 + cb2 * theta + cb3 * pow(theta, 3.0/2.0);
    double ac = ac1 + ac2 * theta + ac3 * pow(theta, 3.0/2.0);
    double bc = bc1 + bc2 * theta + bc3 * pow(theta, 3.0/2.0);
    double cc = cc1 + cc2 * theta + cc3 * pow(theta, 3.0/2.0);
    double ad = ad1 + ad2 * theta + ad3 * pow(theta, 3.0/2.0);
    double bd = bd1 + bd2 * theta + bd3 * pow(theta, 3.0/2.0);
    double cd = cd1 + cd2 * theta + cd3 * pow(theta, 3.0/2.0);
    double xm = Ax + Bx * theta + Cx * pow(theta, 2.0); 

    double g0a = (g0a0 + g0aa1 * theta)/(1 + g0ba1 * theta + g0ba2 * pow(theta, 3.0));
    double g0b = (g0b0 + g0ab1 * sqrt(theta))/(1 + g0bb1 * theta + g0bb2 * pow(theta, 2.0));
    double g0c = (g0c0 + g0ac1 * sqrt(theta) + g0ac2 * pow(theta, 3.0/2.0))/(1 + g0bc1 * theta + g0bc2 * pow(theta, 2.0));
    double g0d = (g0d0 + g0ad1 * sqrt(theta))/(1 + g0bd1 * theta + g0bd2 * pow(theta, 2.0));

    double g00 = 1/2 * (1 + g0a * sqrt(rs) + g0b * rs)/(1 + g0c * rs + g0d * pow(rs, 3.0));

    double a = (aa + ba * rs)/(1 + ca * rs);
    double b = (ab + bb * rs)/(1 + cb * rs);
    double c = (ac + bc * rs)/(1 + cc * rs);
    double d = (ad + bd * rs)/(1 + cd * rs);

    // Assign dx for the free energy derivatives
    double dx = 1e-3; 

    // Assigning derivatives to variables
    double dfxc_rs = derivative_wrt_rs(theta, rs, dx);
    double dfxc_t = derivative_wrt_theta(theta, rs, dx);
    double dfxc_rs2 = second_derivative_wrt_rs(theta, rs, dx);
    double dfxc_t2 = second_derivative_wrt_theta(theta, rs, dx);
    double dfxc_t_rs = mixed_derivative(theta, rs, dx);
    
    const int nx = wvg.size();
    slfc.resize(nx);
    // Loop over the wave vector grid size
    for (int i=0; i<nx; ++i) {
                
        double AF = 1.0/2.0 * (1 + tanh(Eta * (wvg[i] - xm)));

        double GCSR = -(M_PI/12) * lambda * rs * pow(wvg[i], 2.0) * (4 * pow(theta, 2.0) * dfxc_t2 + pow(rs, 2.0) * dfxc_rs2 + 4 * theta * rs * dfxc_t_rs - 2 * theta * dfxc_t - 2 * rs * dfxc_rs);
                
        double Gnnfit = (1 + a * wvg[i] + b * pow(wvg[i], 1.0/2.0))/(1 + c * wvg[i] + d * pow(wvg[i], 1.25) + GCSR);
        
        double Gesa = GCSR * Gnnfit * (1 - AF) + (1 - g00) * AF;
        
        slfc[i] = Gesa;
    }
    // Assign the old slfc values to the slfcOld vector
    slfcOld = slfc;

}

// QMC free energy function constants
void fi(double theta, double &fa, double &fb, double &fd, double &fe, double &fc) {
    // Constant for unit conversion
    const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);

    const double omega = 1;
    const double fb1 = 0.3436902;
    const double fb2 = 7.82159531356;
    const double fb3 = 0.300483986662;
    const double fb4 = 15.8443467125;
    const double fb5 = fb3*pow(3.0/2.0, 1.0/2.0)*omega/lambda;
    const double fc1 = 0.8759442;
    const double fc2 = -0.230130843551;
    const double fd1 = 0.72700876;
    const double fd2 = 2.38264734144;
    const double fd3 = 0.30221237251;
    const double fd4 = 4.39347718395;
    const double fd5 = 0.729951339845;
    const double fe1 = 0.25388214;
    const double fe2 = 0.815795138599;
    const double fe3 = 0.0646844410481;
    const double fe4 = 15.0984620477;
    const double fe5 = 0.230761357474;

    fa = 0.610887 * tanh(pow(theta, -1)) * (0.75 + 3.04363 * pow(theta, 2) - 0.09227 * pow(theta, 3) + 1.7035 * pow(theta, 4)) / (1 + 8.31051 * pow(theta, 2) + 5.1105 * pow(theta, 4));
    fb = tanh(1 / sqrt(theta)) * (fb1 + fb2 * pow(theta, 2) + fb3 * pow(theta, 4)) / (1 + fb4 * pow(theta, 2) + fb5 * pow(theta, 4));
    fd = tanh(1 / sqrt(theta)) * (fd1 + fd2 * pow(theta, 2) + fd3 * pow(theta, 4)) / (1 + fd4 * pow(theta, 2) + fd5 * pow(theta, 4));
    fe = tanh(1 / theta) * (fe1 + fe2 * pow(theta, 2) + fe3 * pow(theta, 4)) / (1 + fe4 * pow(theta, 2) + fe5 * pow(theta, 4));
    fc = (fc1 + fc2 * exp(-(1 / theta))) * fe;
    }

// QMC free energy fitting function
double fxc(double theta, double rs) {
    const double omega = 1;
    double fa, fb, fd, fe, fc;
    fi(theta, fa, fb, fd, fe, fc);
    return -(1 / rs) * (omega * fa + fb * sqrt(rs) + fc * rs) / (1 + fd * sqrt(rs) + fe * rs);
}

// Computation of the derivatives of the free energy functional
double ESA::derivative_wrt_rs(double theta, double rs, double dx) const{
    return (fxc(theta, rs + dx) - fxc(theta, rs - dx)) / (2 * dx);
}

double ESA::derivative_wrt_theta(double theta, double rs, double dx) const{
    return (fxc(theta + dx, rs) - fxc(theta - dx, rs)) / (2 * dx);
}
double ESA::second_derivative_wrt_rs(double theta, double rs, double dx) const{
    return (fxc(theta, rs + dx) - 2 * fxc(theta, rs) + fxc(theta, rs - dx)) / (dx * dx);
}

double ESA::second_derivative_wrt_theta(double theta, double rs, double dx) const{
    return (fxc(theta + dx, rs) - 2 * fxc(theta, rs) + fxc(theta - dx, rs)) / (dx * dx);
}

double ESA::mixed_derivative(double theta, double rs, double dx) const{
    return (fxc(theta + dx, rs + dx) - fxc(theta + dx, rs - dx) - fxc(theta - dx, rs + dx) + fxc(theta - dx, rs - dx)) / (4 * dx * dx);
}