#include "recstls.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "stls.hpp"
#include <cmath>

using namespace std;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;

RecStls::RecStls(const shared_ptr<const RecStlsInput> &in_, const bool verbose_)
    : Rpa(in_, verbose_) {
  ssfInput.resize(wvg.size(), 1.0);
}

const RecStlsInput &RecStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, RecStlsInput>(inPtr);
}

void RecStls::computeStructuralProperties() {
  if (in().getDimension() != dimensionsUtil::Dimension::D3) {
    throwError("REC-STLS is currently implemented only for 3D systems");
  }

  const auto rdfGrid = in().getRdfGrid();
  const auto rdf = in().getRdf();
  if (rdfGrid.empty() || rdf.empty()) {
    throwError("REC-STLS requires a reconstructed RDF grid and values");
  }

  rdfi = make_shared<Interpolator1D>(rdfGrid, rdf);
  if (!rdfi->isValid()) {
    throwError("Unable to interpolate the reconstructed RDF");
  }

  print("Computing reconstructed input structure factor: ");
  computeInputSsf();
  println("Done");

  print("Computing one-shot STLS local field correction: ");
  computeLfc();
  println("Done");

  // Reuse the RPA/STLS structure-factor update, which includes the Hartree-Fock
  // acceleration trick already implemented in Rpa::computeSsf().
  print("Computing output static structure factor: ");
  Rpa::computeSsf();
  println("Done");
}

void RecStls::computeInputSsf() {
  const auto rdfGrid = in().getRdfGrid();
  const double uMin = rdfGrid.front();
  const double uMax = rdfGrid.back();

  for (size_t i = 0; i < wvg.size(); ++i) {
    RecStlsUtil::InputSsf ssfTmp(wvg[i], uMin, uMax, rdfi, itg);
    ssfInput[i] = ssfTmp.get();
  }
}

void RecStls::computeLfc() {
  const int nx = wvg.size();
  const shared_ptr<Interpolator1D> ssfItp =
      make_shared<Interpolator1D>(wvg, ssfInput);
  for (int i = 0; i < nx; ++i) {
    StlsUtil::Slfc lfcTmp(wvg[i], wvg.front(), wvg.back(), ssfItp, itg, inPtr);
    lfc(i, 0) = lfcTmp.get();
  }
}

double RecStlsUtil::InputSsf::get() const {
  auto func = [&](const double &u) -> double { return integrand(u); };
  itg->compute(func, ItgParam(uMin, uMax));
  return 1.0 + 4.0 / (3.0 * M_PI) * itg->getSolution();
}

double RecStlsUtil::InputSsf::integrand(const double &u) const {
  const double u2 = u * u;
  const double gu = rdfi->eval(u) - 1.0;
  if (x == 0.0 || u == 0.0) {
    return u2 * gu;
  }
  return u2 * gu * sin(x * u) / (x * u);
}
