#include "python_wrappers.hpp"
#include "python_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// PyRpaInput
// -----------------------------------------------------------------

bn::ndarray PyRpaInput::getChemicalPotentialGuess(RpaInput &in) {
  return pythonUtil::toNdArray(in.getChemicalPotentialGuess());
}

void PyRpaInput::setChemicalPotentialGuess(RpaInput &in,
                                           const bp::list &muGuess) {
  in.setChemicalPotentialGuess(pythonUtil::toVector(muGuess));
}

// -----------------------------------------------------------------
// PySlfcGuess
// -----------------------------------------------------------------

bn::ndarray PySlfcGuess::getWvg(const StlsInput::SlfcGuess &guess) {
  return pythonUtil::toNdArray(guess.wvg);
}

bn::ndarray PySlfcGuess::getSlfc(const StlsInput::SlfcGuess &guess) {
  return pythonUtil::toNdArray(guess.slfc);
}

void PySlfcGuess::setWvg(StlsInput::SlfcGuess &guess, const bn::ndarray &wvg) {
  guess.wvg = pythonUtil::toVector(wvg);
}

void PySlfcGuess::setSlfc(StlsInput::SlfcGuess &guess,
                          const bn::ndarray &slfc) {
  guess.slfc = pythonUtil::toVector(slfc);
}

// -----------------------------------------------------------------
// PyVSStlsInput
// -----------------------------------------------------------------

bn::ndarray PyVSStlsInput::getAlphaGuess(VSStlsInput &in) {
  return pythonUtil::toNdArray(in.getAlphaGuess());
}

void PyVSStlsInput::setAlphaGuess(VSStlsInput &in, const bp::list &alphaGuess) {
  in.setAlphaGuess(pythonUtil::toVector(alphaGuess));
}

// -----------------------------------------------------------------
// PyQVSStlsInput
// -----------------------------------------------------------------

bn::ndarray PyQVSStlsInput::getAlphaGuess(QVSStlsInput &in) {
  return pythonUtil::toNdArray(in.getAlphaGuess());
}

void PyQVSStlsInput::setAlphaGuess(QVSStlsInput &in,
                                   const bp::list &alphaGuess) {
  in.setAlphaGuess(pythonUtil::toVector(alphaGuess));
}

// -----------------------------------------------------------------
// PyFreeEnergyIntegrand
// -----------------------------------------------------------------

bn::ndarray
PyFreeEnergyIntegrand::getGrid(const VSStlsInput::FreeEnergyIntegrand &fxc) {
  return pythonUtil::toNdArray(fxc.grid);
}

bn::ndarray PyFreeEnergyIntegrand::getIntegrand(
    const VSStlsInput::FreeEnergyIntegrand &fxc) {
  return pythonUtil::toNdArray2D(fxc.integrand);
}

bn::ndarray
PyFreeEnergyIntegrand::getAlpha(const VSStlsInput::FreeEnergyIntegrand &fxc) {
  return pythonUtil::toNdArray(fxc.alpha);
}

void PyFreeEnergyIntegrand::setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
                                    const bn::ndarray &grid) {
  fxc.grid = pythonUtil::toVector(grid);
}

void PyFreeEnergyIntegrand::setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                                         const bn::ndarray &integrand) {
  fxc.integrand = pythonUtil::toDoubleVector(integrand);
}

void PyFreeEnergyIntegrand::setAlpha(VSStlsInput::FreeEnergyIntegrand &fxc,
                                     const bn::ndarray &alpha) {
  fxc.alpha = pythonUtil::toVector(alpha);
}

// -----------------------------------------------------------------
// PyQstlsGuess
// -----------------------------------------------------------------

bn::ndarray PyQstlsGuess::getWvg(const QstlsInput::QstlsGuess &guess) {
  return pythonUtil::toNdArray(guess.wvg);
}

bn::ndarray PyQstlsGuess::getSsf(const QstlsInput::QstlsGuess &guess) {
  return pythonUtil::toNdArray(guess.ssf);
}

bn::ndarray PyQstlsGuess::getAdr(const QstlsInput::QstlsGuess &guess) {
  return pythonUtil::toNdArray2D(guess.adr);
}

int PyQstlsGuess::getMatsubara(const QstlsInput::QstlsGuess &guess) {
  return guess.matsubara;
}

void PyQstlsGuess::setWvg(QstlsInput::QstlsGuess &guess,
                          const bn::ndarray &wvg) {
  guess.wvg = pythonUtil::toVector(wvg);
}

void PyQstlsGuess::setSsf(QstlsInput::QstlsGuess &guess,
                          const bn::ndarray &ssf) {
  guess.ssf = pythonUtil::toVector(ssf);
}

void PyQstlsGuess::setAdr(QstlsInput::QstlsGuess &guess,
                          const bn::ndarray &adr) {
  if (adr.shape(0) == 0) { return; }
  guess.adr = pythonUtil::toVector2D(adr);
}

void PyQstlsGuess::setMatsubara(QstlsInput::QstlsGuess &guess,
                                const int matsubara) {
  guess.matsubara = matsubara;
}

// -----------------------------------------------------------------
// PyRpa
// -----------------------------------------------------------------

int PyRpa::compute(Rpa &rpa) { return rpa.compute(); }

bn::ndarray PyRpa::getIdr(const Rpa &rpa) {
  return pythonUtil::toNdArray2D(rpa.getIdr());
}

bn::ndarray PyRpa::getRdf(const Rpa &rpa, const bn::ndarray &r) {
  return pythonUtil::toNdArray(rpa.getRdf(pythonUtil::toVector(r)));
}

bn::ndarray PyRpa::getSdr(const Rpa &rpa) {
  return pythonUtil::toNdArray(rpa.getSdr());
}

bn::ndarray PyRpa::getSlfc(const Rpa &rpa) {
  return pythonUtil::toNdArray(rpa.getSlfc());
}

bn::ndarray PyRpa::getSsf(const Rpa &rpa) {
  return pythonUtil::toNdArray(rpa.getSsf());
}

bn::ndarray PyRpa::getSsfHF(const Rpa &rpa) {
  return pythonUtil::toNdArray(rpa.getSsfHF());
}

bn::ndarray PyRpa::getWvg(const Rpa &rpa) {
  return pythonUtil::toNdArray(rpa.getWvg());
}

double PyRpa::getUInt(const Rpa &rpa) { return rpa.getUInt(); }

string PyRpa::getRecoveryFileName(const Rpa &rpa) {
  return rpa.getRecoveryFileName();
}

// -----------------------------------------------------------------
// PyStls
// -----------------------------------------------------------------

int PyStls::compute(Stls &stls) { return stls.compute(); }

double PyStls::getError(const Stls &stls) { return stls.getError(); }

bn::ndarray PyStls::getBf(const Stls &stls) {
  return pythonUtil::toNdArray(stls.getBf());
}

// -----------------------------------------------------------------
// PyVSStls
// -----------------------------------------------------------------

int PyVSStls::compute(VSStls &vsstls) { return vsstls.compute(); }

double PyVSStls::getError(const VSStls &vsstls) {
  // NOTE: This is just a place-holder, getError is not yet implemented in
  // VSStls
  if (vsstls.getFreeEnergyIntegrand().empty()) { return -1; }
  return -1;
}

bn::ndarray PyVSStls::getAlpha(const VSStls &vsstls) {
  return pythonUtil::toNdArray(vsstls.getAlpha());
}

bn::ndarray PyVSStls::getFreeEnergyIntegrand(const VSStls &vsstls) {
  return pythonUtil::toNdArray2D(vsstls.getFreeEnergyIntegrand());
}

bn::ndarray PyVSStls::getFreeEnergyGrid(const VSStls &vsstls) {
  return pythonUtil::toNdArray(vsstls.getFreeEnergyGrid());
}

// -----------------------------------------------------------------
// PyQstls
// -----------------------------------------------------------------

int PyQstls::compute(Qstls &qstls) { return qstls.compute(); }

double PyQstls::getError(const Qstls &qstls) { return qstls.getError(); }

bn::ndarray PyQstls::getAdr(const Qstls &qstls) {
  return pythonUtil::toNdArray2D(qstls.getAdr());
}

// -----------------------------------------------------------------
// PyQVSStls
// -----------------------------------------------------------------

int PyQVSStls::compute(QVSStls &qvsstls) { return qvsstls.compute(); }

double PyQVSStls::getError(const QVSStls &qvsstls) {
  // NOTE: This is just a place-holder, getError is not yet implemented in
  // QVSStls
  if (qvsstls.getFreeEnergyIntegrand().empty()) { return -1; }
  return -1;
}

bn::ndarray PyQVSStls::getAlpha(const QVSStls &qvsstls) {
  return pythonUtil::toNdArray(qvsstls.getAlpha());
}

bn::ndarray PyQVSStls::getAdr(const QVSStls &qvsstls) {
  return pythonUtil::toNdArray2D(qvsstls.getAdr());
}

bn::ndarray PyQVSStls::getFreeEnergyIntegrand(const QVSStls &qvsstls) {
  return pythonUtil::toNdArray2D(qvsstls.getFreeEnergyIntegrand());
}

bn::ndarray PyQVSStls::getFreeEnergyGrid(const QVSStls &qvsstls) {
  return pythonUtil::toNdArray(qvsstls.getFreeEnergyGrid());
}

// -----------------------------------------------------------------
// PyThermo
// -----------------------------------------------------------------

bn::ndarray PyThermo::computeRdf(const bn::ndarray &rIn,
                                 const bn::ndarray &wvgIn,
                                 const bn::ndarray &ssfIn) {
  const vector<double> &r = pythonUtil::toVector(rIn);
  const vector<double> &wvg = pythonUtil::toVector(wvgIn);
  const vector<double> &ssf = pythonUtil::toVector(ssfIn);
  return pythonUtil::toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
}

double PyThermo::computeInternalEnergy(const bn::ndarray &wvgIn,
                                       const bn::ndarray &ssfIn,
                                       const double &coupling) {
  const vector<double> &wvg = pythonUtil::toVector(wvgIn);
  const vector<double> &ssf = pythonUtil::toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
}

double PyThermo::computeFreeEnergy(const bn::ndarray &gridIn,
                                   const bn::ndarray &rsuIn,
                                   const double &coupling) {
  const vector<double> &grid = pythonUtil::toVector(gridIn);
  const vector<double> &rsu = pythonUtil::toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}
