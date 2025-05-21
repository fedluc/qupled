#include "python_wrappers.hpp"
#include "python_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// PyInput
// -----------------------------------------------------------------

bn::ndarray PyInput::getChemicalPotentialGuess(Input &in) {
  return pythonUtil::toNdArray(in.getChemicalPotentialGuess());
}

void PyInput::setChemicalPotentialGuess(Input &in, const bp::list &muGuess) {
  in.setChemicalPotentialGuess(pythonUtil::toVector(muGuess));
}

// -----------------------------------------------------------------
// PyGuess
// -----------------------------------------------------------------

bn::ndarray PyGuess::getWvg(const Guess &guess) {
  return pythonUtil::toNdArray(guess.wvg);
}

bn::ndarray PyGuess::getSsf(const Guess &guess) {
  return pythonUtil::toNdArray(guess.ssf);
}

bn::ndarray PyGuess::getLfc(const Guess &guess) {
  return pythonUtil::toNdArray2D(guess.lfc);
}

void PyGuess::setWvg(Guess &guess, const bn::ndarray &wvg) {
  guess.wvg = pythonUtil::toVector(wvg);
}

void PyGuess::setSsf(Guess &guess, const bn::ndarray &ssf) {
  guess.ssf = pythonUtil::toVector(ssf);
}

void PyGuess::setLfc(Guess &guess, const bn::ndarray &lfc) {
  if (lfc.shape(0) == 0) { return; }
  guess.lfc = pythonUtil::toVector2D(lfc);
}

// -----------------------------------------------------------------
// PyVSInput
// -----------------------------------------------------------------

bn::ndarray PyVSInput::getAlphaGuess(VSInput &in) {
  return pythonUtil::toNdArray(in.getAlphaGuess());
}

void PyVSInput::setAlphaGuess(VSInput &in, const bp::list &alphaGuess) {
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

void PyFreeEnergyIntegrand::setGrid(VSStlsInput::FreeEnergyIntegrand &fxc,
                                    const bn::ndarray &grid) {
  fxc.grid = pythonUtil::toVector(grid);
}

void PyFreeEnergyIntegrand::setIntegrand(VSStlsInput::FreeEnergyIntegrand &fxc,
                                         const bn::ndarray &integrand) {
  fxc.integrand = pythonUtil::toDoubleVector(integrand);
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

// -----------------------------------------------------------------
// PyDatabaseInfo
// -----------------------------------------------------------------

string PyDatabaseInfo::getName(const DatabaseInfo &db) { return db.name; }

string PyDatabaseInfo::getRunTableName(const DatabaseInfo &db) {
  return db.runTableName;
}

int PyDatabaseInfo::getRunId(const DatabaseInfo &db) { return db.runId; }

void PyDatabaseInfo::setName(DatabaseInfo &db, const string &name) {
  db.name = name;
}

void PyDatabaseInfo::setRunTableName(DatabaseInfo &db,
                                     const string &runTableName) {
  db.runTableName = runTableName;
}

void PyDatabaseInfo::setRunId(DatabaseInfo &db, const int runId) {
  db.runId = runId;
}