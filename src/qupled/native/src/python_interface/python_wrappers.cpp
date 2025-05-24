#include "python_wrappers.hpp"
#include "python_util.hpp"

using namespace std;

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