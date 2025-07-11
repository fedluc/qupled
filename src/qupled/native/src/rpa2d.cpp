// #include "rpa.hpp"
// #include "chemical_potential.hpp"
// #include "input.hpp"
// #include "mpi_util.hpp"
// #include "numerics.hpp"
// #include "thermo_util.hpp"
// #include <cmath>

// using namespace std;
// using namespace thermoUtil;
// using namespace MPIUtil;
// using ItgParam = Integrator1D::Param;
// using ItgType = Integrator1D::Type;

// // Compute 2D static structure factor at finite temperature
// void Rpa::computeSsfFinite2D() {
//   const double Theta = in().getDegeneracy();
//   const double rs = in().getCoupling();
//   const size_t nx = wvg.size();
//   for (size_t i = 0; i < nx; ++i) {
//     RpaUtil::Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], lfc[i], idr[i]);
//     ssf[i] = ssfTmp.get();
//   }
// }

// // Compute 2D static local field correction
// void Rpa::computeLfc2D() {
//   assert(lfc.size() == wvg.size());
//   for (auto &s : lfc) {
//     s = 0;
//   }
// }

// // Get at finite temperature for any 2D scheme
// double RpaUtil::Ssf::get() const {
//   if (rs == 0.0) return ssfHF;
//   if (x == 0.0) return 0.0;
//   const double isStatic = lfc.size() == 1;
//   double suml = 0.0;
//   for (size_t l = 0; l < idr.size(); ++l) {
//     const double &idrl = idr[l];
//     const double &lfcl = (isStatic) ? lfc[0] : lfc[l];
//     const double denom = 1.0 + ip * idrl * (1 - lfcl);
//     const double f = idrl * idrl * (1 - lfcl) / denom;
//     suml += (l == 0) ? f : 2 * f;
//   }
//   return ssfHF - 1.5 * ip * Theta * suml;
// }