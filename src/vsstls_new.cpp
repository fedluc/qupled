#include "vsstls_new.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

double VSStlsNew::computeAlpha() {
  // Compute the free energy integrand
  thermoProp.compute();
  // Free energy
  const vector<double> freeEnergyData = thermoProp.getFreeEnergyData();
  const double &fxc = freeEnergyData[0];
  const double &fxcr = freeEnergyData[1];
  const double &fxcrr = freeEnergyData[2];
  const double &fxct = freeEnergyData[3];
  const double &fxctt = freeEnergyData[4];
  const double &fxcrt = freeEnergyData[5];
  // Internal energy
  const vector<double> internalEnergyData = thermoProp.getInternalEnergyData();
  const double &uint = internalEnergyData[0];
  const double &uintr = internalEnergyData[1];
  const double &uintt = internalEnergyData[2];
  // Alpha
  double numer = 2 * fxc - (1.0 / 6.0) * fxcrr + (4.0 / 3.0) * fxcr;
  double denom = uint + (1.0 / 3.0) * uintr;
  if (in.getDegeneracy() > 0.0) {
    numer += -(2.0 / 3.0) * fxctt - (2.0 / 3.0) * fxcrt + (1.0 / 3.0) * fxct;
    denom += (2.0 / 3.0) * uintt;
  }
  return numer / denom;
}

void VSStlsNew::updateSolution() {
  // Update the structural properties used for output
  slfc = thermoProp.getSlfc();
  ssf = thermoProp.getSsf();
}

void VSStlsNew::initScheme() {
  Rpa::init();
}

void VSStlsNew::initFreeEnergyIntegrand() {
  if (!thermoProp.isFreeEnergyIntegrandIncomplete()) { return; }
  if (verbose) {
    printf("Missing points in the free energy integrand: subcalls will be "
           "performed to collect the necessary data\n");
  }
  if (verbose) {
    printf("-----------------------------------------------------------------"
           "----------\n");
  }
  VSStlsInput inTmp = in;
  while (thermoProp.isFreeEnergyIntegrandIncomplete()) {
    const double rs = thermoProp.getFirstUnsolvedStatePoint();
    if (verbose) { printf("Subcall: solving VS scheme for rs = %.5f:\n", rs); }
    inTmp.setCoupling(rs);
    VSStlsNew scheme(inTmp, thermoProp);
    scheme.compute();
    thermoProp.copyFreeEnergyIntegrand(scheme.getThermoProp());
    if (verbose) {
      printf("Done\n");
      printf("-----------------------------------------------------------------"
             "----------\n");
    }
  }
}

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

void StructPropNew::doIterations() {
 //  const auto &in = csr[0].getInput();
//   const int maxIter = in.getNIter();
//   const int ompThreads = in.getNThreads();
//   const double minErr = in.getErrMin();
//   double err = 1.0;
//   int counter = 0;
//   // Define initial guess
//   for (auto &c : csr) {
//     c.initialGuess();
//   }
//   // Iteration to solve for the structural properties
//   const bool useOMP = ompThreads > 1;
//   while (counter < maxIter + 1 && err > minErr) {
// // Compute new solution and error
// #pragma omp parallel num_threads(ompThreads) if (useOMP)
//     {
// #pragma omp for
//       for (auto &c : csr) {
//         c.computeSsf();
//         c.computeSlfcStls();
//       }
// #pragma omp for
//       for (size_t i = 0; i < csr.size(); ++i) {
//         auto &c = csr[i];
//         c.computeSlfc();
//         if (i == RS_THETA) { err = c.computeError(); }
//         c.updateSolution();
//       }
//     }
//     counter++;
//   }
//   if (verbose) {
//     printf("Alpha = %.5e, Residual error "
//            "(structural properties) = %.5e\n",
//            csr[RS_THETA].getAlpha(),
//            err);
//   }
}
