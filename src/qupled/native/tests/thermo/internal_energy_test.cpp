#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "thermo/internal_energy.hpp"
#include "util/numerics.hpp"

TEST(InternalEnergyTest, IsZeroForUnitStructureFactorIn2D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const InternalEnergy u2d(
      2.0, wvg.front(), wvg.back(), itp, itg, dimensionsUtil::Dimension::D2);
  EXPECT_NEAR(u2d.get(), 0.0, 1.0e-10);
}

TEST(InternalEnergyTest, IsZeroForUnitStructureFactorIn3D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const InternalEnergy u3d(
      2.0, wvg.front(), wvg.back(), itp, itg, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(u3d.get(), 0.0, 1.0e-10);
}

TEST(InternalEnergyTest, MatchesClosedFormForConstantOffsetStructureFactor) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 2.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const InternalEnergy u2d(
      2.0, 0.0, 2.0, itp, itg, dimensionsUtil::Dimension::D2);
  EXPECT_NEAR(u2d.get(), std::sqrt(2.0) / 2.0, 1.0e-6);

  const InternalEnergy u3d(
      2.0, 0.0, 2.0, itp, itg, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(u3d.get(), 1.0 / (M_PI * numUtil::lambda), 1.0e-6);
}
