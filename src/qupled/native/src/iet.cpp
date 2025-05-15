#include "iet.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include <format>

using namespace std;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

// Initialize basic properties
template <typename IetInput>
void Iet<IetInput>::init() {
  print("Computing bridge function adder: ");
  computeBf();
  println("Done");
}

// Compute bridge function
template <typename IetInput>
void Iet<IetInput>::computeBf() {
  const size_t nx = bf.size();
  const shared_ptr<Integrator1D> itgF =
      make_shared<Integrator1D>(ItgType::FOURIER, 1e-10);
  assert(bf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    IetUtil::BridgeFunction bfTmp(in->getTheory(),
                                  in->getMapping(),
                                  in->getCoupling(),
                                  in->getDegeneracy(),
                                  wvg[i],
                                  itgF);
    bf[i] = bfTmp.get();
  }
}

// Read initial guess from input
template <typename IetInput>
bool Iet<IetInput>::initialGuessFromInput(Vector2D &lfc) {
  const auto &guess = in->getGuess();
  const int nx = lfc.size(0);
  const int nl = lfc.size(1);
  const int nx_ = guess.lfc.size(0);
  const int nl_ = guess.lfc.size(1);
  const double xMax = (guess.wvg.empty()) ? 0.0 : guess.wvg.back();
  vector<Interpolator1D> itp(nl_);
  for (int l = 0; l < nl_; ++l) {
    vector<double> tmp(nx_);
    for (int i = 0; i < nx_; ++i) {
      tmp[i] = guess.lfc(i, l);
    }
    itp[l].reset(guess.wvg[0], tmp[0], nx_);
    if (!itp[l].isValid()) { return false; }
  }
  for (int i = 0; i < nx; ++i) {
    const double &x = wvg[i];
    if (x > xMax) {
      lfc.fill(i, 0.0);
      continue;
    }
    for (int l = 0; l < nl; ++l) {
      lfc(i, l) = (l < nl_) ? itp[l].eval(x) : 0.0;
    }
  }
  return true;
}

// Explicit instantiation
template class Iet<StlsIetInput>;
template class Iet<QstlsIetInput>;

// -----------------------------------------------------------------
// BridgeFunction class
// -----------------------------------------------------------------

double IetUtil::BridgeFunction::get() const {
  if (theory == "STLS-HNC" || theory == "QSTLS-HNC") { return hnc(); }
  if (theory == "STLS-IOI" || theory == "QSTLS-IOI") { return ioi(); }
  if (theory == "STLS-LCT" || theory == "QSTLS-LCT") { return lct(); }
  throwError("Unknown theory to compute the bridge function term");
  return numUtil::Inf;
}

double IetUtil::BridgeFunction::couplingParameter() const {
  const double fact = 2 * lambda * lambda * rs;
  if (mapping == "sqrt") { return fact / sqrt(1 + Theta * Theta); }
  if (mapping == "linear") { return fact / (1 + Theta); }
  if (Theta != 0.0) { return fact / Theta; }
  throwError("The standard iet mapping cannot be used in the "
             "ground state");
  return numUtil::Inf;
}

double IetUtil::BridgeFunction::hnc() const { return 0.0; }

double IetUtil::BridgeFunction::ioi() const {
  const double l2 = lambda * lambda;
  const double l3 = l2 * lambda;
  const double l4 = l3 * lambda;
  const double l5 = l4 * lambda;
  const double l6 = l5 * lambda;
  const double l7 = l6 * lambda;
  const double l8 = l7 * lambda;
  const double Gamma = couplingParameter();
  const double lnG = log(Gamma);
  const double lnG2 = lnG * lnG;
  const double b0 = 0.258 - 0.0612 * lnG + 0.0123 * lnG2 - 1.0 / Gamma;
  const double b1 = 0.0269 + 0.0318 * lnG + 0.00814 * lnG2;
  if (b0 / b1 <= 0.0 || Gamma < 5.25 || Gamma > 171.8) {
    const string msg = format("The IET schemes cannot be applied "
                              "to this state point because Gamma = {:.8f} "
                              "falls outside the range of validty of the "
                              "bridge function parameterization\n",
                              Gamma);
    throwError(msg);
  }
  const double c1 = 0.498 - 0.280 * lnG + 0.0294 * lnG2;
  const double c2 = -0.412 + 0.219 * lnG - 0.0251 * lnG2;
  const double c3 = 0.0988 - 0.0534 * lnG + 0.00682 * lnG2;
  const double b02 = b0 * b0;
  const double b03 = b02 * b0;
  const double b04 = b03 * b0;
  const double b05 = b04 * b0;
  const double b06 = b05 * b0;
  const double b07 = b06 * b0;
  const double b08 = b07 * b0;
  const double b12 = b1 * b1;
  const double b13 = b12 * b1;
  const double b14 = b13 * b1;
  const double b15 = b14 * b1;
  const double b16 = b15 * b1;
  const double b17 = b16 * b1;
  const double b18 = b17 * b1;
  const double b02_b12 = b02 / b12;
  const double b03_b13 = b03 / b13;
  const double b04_b14 = b04 / b14;
  const double b05_b15 = b05 / b15;
  const double b06_b16 = b06 / b16;
  const double b07_b17 = b07 / b17;
  const double b08_b18 = b08 / b18;
  const double fact = sqrt(M_PI) / (4.0 * l2) * pow(b0 / b1, 1.5);
  const double q2 = x * x;
  const double q3 = q2 * x;
  const double q4 = q3 * x;
  const double q5 = q4 * x;
  const double q6 = q5 * x;
  const double q7 = q6 * x;
  const double q8 = q7 * x;
  const double bf1 =
      -b0
      + c1 / 16.0
            * (60.0 * b02_b12 - 20.0 * b03_b13 * q2 / l2 + b04_b14 * q4 / l4);
  const double bf2 = c2 / 64.0
                     * (840.0 * b03_b13 - 420.0 * b04_b14 * q2 / l2
                        + 42.0 * b05_b15 * q4 / l4 - b06_b16 * q6 / l6);
  ;
  const double bf3 = c3 / 256.0
                     * (15120.0 * b04_b14 - 10080.0 * b05_b15 * q2 / l2
                        + 1512.0 * b06_b16 * q4 / l4 - 72.0 * b07_b17 * q6 / l6
                        + b08_b18 * q8 / l8);
  return fact * q2 * (bf1 + bf2 + bf3) * exp(-b0 * q2 / (4.0 * b1 * l2));
}

double IetUtil::BridgeFunction::lct() const {
  const double Gamma = couplingParameter();
  auto func = [&](const double &r) -> double { return lctIntegrand(r, Gamma); };
  itg->compute(func, ItgParam(x / lambda));
  return itg->getSolution() * (x / lambda) / Gamma;
  return 0.0;
}

double IetUtil::BridgeFunction::lctIntegrand(const double &r,
                                             const double &Gamma) const {
  if (Gamma < 5.0) {
    const string msg = format("The IET schemes cannot be applied "
                              "to this state point because Gamma = {:.3f} "
                              "falls outside the range of validty of the "
                              "bridge function parameterization\n",
                              Gamma);
    throwError(msg);
  }
  const double Gamma1_6 = pow(Gamma, 1. / 6.);
  const double lnG = log(Gamma);
  const double lnG2 = lnG * lnG;
  const double lnG3 = lnG2 * lnG;
  const double lnG4 = lnG3 * lnG;
  const double a0 =
      Gamma * (0.076912 - 0.10465 * lnG + 0.0056629 * lnG2 + 0.00025656 * lnG3);
  const double a2 =
      Gamma * (0.068045 - 0.036952 * lnG + 0.048818 * lnG2 - 0.0048985 * lnG3);
  const double a3 =
      Gamma * (-0.30231 + 0.30457 * lnG - 0.11424 * lnG2 + 0.0095993 * lnG3);
  const double a4 =
      Gamma * (0.25111 - 0.26800 * lnG + 0.082268 * lnG2 - 0.0064960 * lnG3);
  const double a5 =
      Gamma * (-0.061894 + 0.066811 * lnG - 0.019140 * lnG2 + 0.0014743 * lnG3);
  const double c0 = Gamma
                    * (0.25264 - 0.31615 * lnG + 0.13135 * lnG2
                       - 0.023044 * lnG3 + 0.0014666 * lnG4);
  const double c1 = Gamma1_6
                    * (-12.665 + 20.802 * lnG - 9.6296 * lnG2 + 1.7889 * lnG3
                       - 0.11810 * lnG4);
  const double c2 = Gamma1_6
                    * (15.285 - 14.076 * lnG + 5.7558 * lnG2 - 1.0188 * lnG3
                       + 0.06551 * lnG4);
  const double c3 = Gamma1_6
                    * (35.330 - 40.727 * lnG + 16.690 * lnG2 - 2.8905 * lnG3
                       + 0.18243 * lnG4);
  const double r2 = r * r;
  const double r3 = r2 * r;
  const double r4 = r3 * r;
  const double r5 = r4 * r;
  const double rshift = r - 1.44;
  const double bsr = a0 + a2 * r2 + a3 * r3 + a4 * r4 + a5 * r5;
  const double blr = c0 * exp(-c1 * rshift) * exp(-0.3 * r2)
                     * (cos(c2 * rshift) + c3 * exp(-3.5 * rshift));
  const double sf = 0.5 * (1.0 + erf(5.0 * (r - 1.50)));
  return r * ((1 - sf) * bsr + sf * blr);
}
