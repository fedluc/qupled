#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "util.hpp"
#include "python_wrappers.hpp"

namespace vp = vecUtil::python;
using namespace thermoUtil;

// -----------------------------------------------------------------
// PyRpa
// -----------------------------------------------------------------

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

