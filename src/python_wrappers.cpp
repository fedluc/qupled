#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "util.hpp"
#include "numerics.hpp"
#include "input.hpp"
#include "python_wrappers.hpp"

namespace vp = vecUtil::python;
using namespace std;
using namespace thermoUtil;

// -----------------------------------------------------------------
// PyRpaInput
// -----------------------------------------------------------------

bn::ndarray PyRpaInput::getChemicalPotentialGuess(RpaInput &in){
  return vp::toNdArray(in.getChemicalPotentialGuess());
}
  
void PyRpaInput::setChemicalPotentialGuess(RpaInput &in,
					const bp::list &muGuess){
  in.setChemicalPotentialGuess(vp::toVector(muGuess));
}

// -----------------------------------------------------------------
// PySlfcGuess
// -----------------------------------------------------------------

bn::ndarray PySlfcGuess::getWvg(const StlsInput::SlfcGuess &guess){
  return vp::toNdArray(guess.wvg);
}

bn::ndarray PySlfcGuess::getSlfc(const StlsInput::SlfcGuess &guess){
  return vp::toNdArray(guess.slfc);
}

void PySlfcGuess::setWvg(StlsInput::SlfcGuess &guess,
			 const bn::ndarray& wvg) {
  guess.wvg = vp::toVector(wvg);
}

void PySlfcGuess::setSlfc(StlsInput::SlfcGuess &guess,
			  const bn::ndarray& slfc) {
  guess.slfc = vp::toVector(slfc);
}

// -----------------------------------------------------------------
// PyVSStlsInput
// -----------------------------------------------------------------

bn::ndarray PyVSStlsInput::getAlphaGuess(VSStlsInput &in){
  return vp::toNdArray(in.getAlphaGuess());
}
  
void PyVSStlsInput::setAlphaGuess(VSStlsInput &in,
				  const bp::list &alphaGuess){
  in.setAlphaGuess(vp::toVector(alphaGuess));
}

// -----------------------------------------------------------------
// PyQVSStlsInput
// -----------------------------------------------------------------

bn::ndarray PyQVSStlsInput::getAlphaGuess(QVSStlsInput &in){
  return vp::toNdArray(in.getAlphaGuess());
}
  
void PyQVSStlsInput::setAlphaGuess(QVSStlsInput &in,
				   const bp::list &alphaGuess){
  in.setAlphaGuess(vp::toVector(alphaGuess));
}


// -----------------------------------------------------------------
// PyFreeEnergyIntegrand
// -----------------------------------------------------------------

bn::ndarray PyFreeEnergyIntegrand::getGrid(const VSStlsInput::FreeEnergyIntegrand &fxc){
  return vp::toNdArray(fxc.grid);
}

bn::ndarray PyFreeEnergyIntegrand::getIntegrand(const VSStlsInput::FreeEnergyIntegrand&fxc){
  return vp::toNdArray2D(fxc.integrand);
}

void PyFreeEnergyIntegrand::setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
				    const bn::ndarray& grid) {
  fxc.grid = vp::toVector(grid);
}

void PyFreeEnergyIntegrand::setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
					 const bn::ndarray& integrand) {
  fxc.integrand = vp::toDoubleVector(integrand);
}

// -----------------------------------------------------------------
// PyQstlsGuess
// -----------------------------------------------------------------

bn::ndarray PyQstlsGuess::getWvg(const QstlsInput::QstlsGuess &guess){
  return vp::toNdArray(guess.wvg);
}

bn::ndarray PyQstlsGuess::getSsf(const QstlsInput::QstlsGuess &guess){
  return vp::toNdArray(guess.ssf);
}

bn::ndarray PyQstlsGuess::getAdr(const QstlsInput::QstlsGuess &guess){
  return vp::toNdArray2D(guess.adr);
}

int PyQstlsGuess::getMatsubara(const QstlsInput::QstlsGuess &guess){
  return guess.matsubara;
}

void PyQstlsGuess::setWvg(QstlsInput::QstlsGuess &guess,
			 const bn::ndarray& wvg) {
  guess.wvg = vp::toVector(wvg);
}

void PyQstlsGuess::setSsf(QstlsInput::QstlsGuess &guess,
			  const bn::ndarray& ssf) {
  guess.ssf = vp::toVector(ssf);
}

void PyQstlsGuess::setAdr(QstlsInput::QstlsGuess &guess,
			  const bn::ndarray& adr) {
  if (adr.shape(0) == 0) {
    return;
  }
  guess.adr = vp::toVector2D(adr);
}

void PyQstlsGuess::setMatsubara(QstlsInput::QstlsGuess &guess,
				const int matsubara){
  guess.matsubara = matsubara;
}

// -----------------------------------------------------------------
// PyRpa
// -----------------------------------------------------------------

int PyRpa::compute(Rpa& rpa) {
  return rpa.compute();
}

bn::ndarray PyRpa::getIdr(const Rpa& rpa) {
  return vp::toNdArray2D(rpa.getIdr());
}

bn::ndarray PyRpa::getRdf(const Rpa& rpa,
			  const bn::ndarray &r) {
  return vp::toNdArray(rpa.getRdf(vp::toVector(r)));
}

bn::ndarray PyRpa::getSdr(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSdr());
}

bn::ndarray PyRpa::getSlfc(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSlfc());
}

bn::ndarray PyRpa::getSsf(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSsf());
}

bn::ndarray PyRpa::getSsfHF(const Rpa& rpa) {
  return vp::toNdArray(rpa.getSsfHF());
}

bn::ndarray PyRpa::getWvg(const Rpa& rpa) {
  return vp::toNdArray(rpa.getWvg());
}

double PyRpa::getUInt(const Rpa& rpa) {
  return rpa.getUInt();
}

string PyRpa::getRecoveryFileName(const Rpa& rpa) {
  return rpa.getRecoveryFileName();
}

// -----------------------------------------------------------------
// PyStls
// -----------------------------------------------------------------

int PyStls::compute(Stls& stls) {
  return stls.compute();
}

bn::ndarray PyStls::getBf(const Stls& stls) {
  return vp::toNdArray(stls.getBf());
}

// -----------------------------------------------------------------
// PyVSStls
// -----------------------------------------------------------------

int PyVSStls::compute(VSStls& vsstls) {
  return vsstls.compute();
}

bn::ndarray PyVSStls::getFreeEnergyIntegrand(const VSStls &vsstls){
  return vp::toNdArray2D(vsstls.getFreeEnergyIntegrand());
}

bn::ndarray PyVSStls::getFreeEnergyGrid(const VSStls &vsstls){
  return vp::toNdArray(vsstls.getFreeEnergyGrid());
}

// -----------------------------------------------------------------
// PyQVSStls
// -----------------------------------------------------------------

int PyQVSStls::compute(QVSStls& qvsstls) {
  return qvsstls.compute();
}

// -----------------------------------------------------------------
// PyQstls
// -----------------------------------------------------------------

int PyQstls::compute(Qstls& qstls) {
  return qstls.compute();
}

bn::ndarray PyQstls::getAdr(const Qstls& qstls) {
  return vp::toNdArray2D(qstls.getAdr());
}

// -----------------------------------------------------------------
// PyThermo
// -----------------------------------------------------------------

bn::ndarray PyThermo::computeRdf(const bn::ndarray &rIn,
				 const bn::ndarray &wvgIn,
				 const bn::ndarray &ssfIn) {
  const vector<double> &r = vp::toVector(rIn);
  const vector<double> &wvg = vp::toVector(wvgIn);
  const vector<double> &ssf = vp::toVector(ssfIn);
  return vp::toNdArray(thermoUtil::computeRdf(r, wvg, ssf));
}

double PyThermo::computeInternalEnergy(const bn::ndarray &wvgIn,
				       const bn::ndarray &ssfIn,
				       const double &coupling) {
  const vector<double> &wvg = vp::toVector(wvgIn);
  const vector<double> &ssf = vp::toVector(ssfIn);
  return thermoUtil::computeInternalEnergy(wvg, ssf, coupling);
}

double PyThermo::computeFreeEnergy(const bn::ndarray &gridIn,
				   const bn::ndarray &rsuIn,
				   const double &coupling) {
  const vector<double> &grid = vp::toVector(gridIn);
  const vector<double> &rsu = vp::toVector(rsuIn);
  return thermoUtil::computeFreeEnergy(grid, rsu, coupling);
}

