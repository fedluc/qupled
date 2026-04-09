#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "schemes/input.hpp"
#include "thermo/chemical_potential.hpp"
#include "thermo/free_energy.hpp"
#include "thermo/internal_energy.hpp"
#include "thermo/rdf.hpp"
#include "util/dimensions_util.hpp"
#include "util/numerics.hpp"

TEST(ThermoTest, FreeEnergyIntegratesInterpolatedCurveAndNormalizationFlag) {
  const std::vector<double> rs_grid{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> rsu{2.0, 2.0, 2.0, 2.0};
  auto itp = std::make_shared<Interpolator1D>(rs_grid, rsu);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const FreeEnergy non_norm(2.0, itp, itg, false);
  EXPECT_NEAR(non_norm.get(), 4.0, 1.0e-6);

  const FreeEnergy norm(2.0, itp, itg, true);
  EXPECT_NEAR(norm.get(), 1.0, 1.0e-6);

  const FreeEnergy zero_norm(0.0, itp, itg, true);
  EXPECT_TRUE(std::isinf(zero_norm.get()));
  EXPECT_LT(zero_norm.get(), 0.0);
}

TEST(ThermoTest, InternalEnergyIsZeroForUnitStaticStructureFactor) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const InternalEnergy u2d(2.0,
                           wvg.front(),
                           wvg.back(),
                           itp,
                           itg,
                           dimensionsUtil::Dimension::D2);
  EXPECT_NEAR(u2d.get(), 0.0, 1.0e-10);

  const InternalEnergy u3d(2.0,
                           wvg.front(),
                           wvg.back(),
                           itp,
                           itg,
                           dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(u3d.get(), 0.0, 1.0e-10);
}

TEST(ThermoTest, RdfStaysAtOneWhenStructureFactorIsOne) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-8);
  auto itgf = std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  Rdf rdf2d(0.3, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D2);
  EXPECT_NEAR(rdf2d.get(), 1.0, 1.0e-6);

  Rdf rdf3d_r0(0.0, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(rdf3d_r0.get(), 1.0, 1.0e-6);

  Rdf rdf3d(1.2, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(rdf3d.get(), 1.0, 1.0e-5);
}

TEST(ThermoTest, ChemicalPotentialMatchesClosedFormAndRootCondition) {
  auto in2d = std::make_shared<Input>();
  in2d->setDegeneracy(1.25);
  ChemicalPotential cp2d(in2d);
  cp2d.compute(dimensionsUtil::Dimension::D2);
  const double expected2d = std::log(std::exp(1.0 / 1.25) - 1.0);
  EXPECT_NEAR(cp2d.get(), expected2d, 1.0e-12);

  auto in3d = std::make_shared<Input>();
  in3d->setDegeneracy(2.0);
  in3d->setChemicalPotentialGuess({-10.0, 10.0});
  ChemicalPotential cp3d(in3d);
  cp3d.compute(dimensionsUtil::Dimension::D3);
  const double residual =
      SpecialFunctions::fermiDirac12(cp3d.get())
      - 2.0 / (3.0 * std::pow(in3d->getDegeneracy(), 1.5));
  EXPECT_NEAR(residual, 0.0, 1.0e-8);
}
