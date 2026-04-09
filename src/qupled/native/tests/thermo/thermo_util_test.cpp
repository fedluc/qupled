#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "fixtures/tolerance.hpp"
#include "thermo/thermo_util.hpp"

TEST(ThermoUtilTest, ComputeInternalEnergyReturnsZeroForUnitStructureFactorIn2D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);

  EXPECT_NEAR(
      thermoUtil::computeInternalEnergy(wvg, ssf, 1.0, dimensionsUtil::Dimension::D2),
      0.0,
      testFixtures::kTightTol);
}

TEST(ThermoUtilTest, ComputeInternalEnergyReturnsZeroForUnitStructureFactorIn3D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);

  EXPECT_NEAR(
      thermoUtil::computeInternalEnergy(wvg, ssf, 1.0, dimensionsUtil::Dimension::D3),
      0.0,
      testFixtures::kTightTol);
}

TEST(ThermoUtilTest, ComputeFreeEnergyOverloadsHandleNormalizationAndRangeChecks) {
  const std::vector<double> grid{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> rsu{2.0, 2.0, 2.0, 2.0};

  EXPECT_NEAR(thermoUtil::computeFreeEnergy(grid, rsu, 2.0, false), 4.0, 1.0e-6);
  EXPECT_NEAR(thermoUtil::computeFreeEnergy(grid, rsu, 2.0, true), 1.0, 1.0e-6);
  EXPECT_NEAR(thermoUtil::computeFreeEnergy(grid, rsu, 2.0), 1.0, 1.0e-6);
  EXPECT_THROW(thermoUtil::computeFreeEnergy(grid, rsu, 3.1), std::runtime_error);
}

TEST(ThermoUtilTest, ComputeRdfReturnsExpectedShapeAndBaselineValues) {
  const std::vector<double> r{0.0, 0.3, 1.0};
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);

  const auto rdf2d = thermoUtil::computeRdf(r, wvg, ssf, dimensionsUtil::Dimension::D2);
  ASSERT_EQ(rdf2d.size(), r.size());
  for (double v : rdf2d) {
    EXPECT_NEAR(v, 1.0, 1.0e-5);
  }
}

TEST(ThermoUtilTest, ComputeItcfReturnsEmptyForEmptyAxes) {
  auto in = testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  const Vector2D idr(1, 1);
  const Vector2D lfc(1, 1);

  EXPECT_TRUE(thermoUtil::computeItcf(in, {}, {0.0}, 0.0, idr, lfc).empty());
  EXPECT_TRUE(thermoUtil::computeItcf(in, {0.0, 1.0}, {}, 0.0, idr, lfc).empty());
}

TEST(ThermoUtilTest, ComputeItcfUsesHartreeFockBypassWhenTheoryIsHf) {
  const std::vector<double> wvg{0.0, 1.0};
  const std::vector<double> tau{0.0};
  const Vector2D idr(2, 1);
  const Vector2D lfcMismatched(1, 1);

  auto hfInput = testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  const Vector2D hfRes = thermoUtil::computeItcf(hfInput, wvg, tau, 0.0, idr, lfcMismatched);
  ASSERT_EQ(hfRes.size(0), 2u);
  ASSERT_EQ(hfRes.size(1), 1u);
  EXPECT_NEAR(hfRes(0, 0), 0.0, testFixtures::kTightTol);
  EXPECT_NEAR(hfRes(1, 0), 11.0 / 16.0, testFixtures::kTightTol);
}

TEST(ThermoUtilTest, ComputeItcfValidatesLfcShapeForInteractingTheory) {
  const std::vector<double> wvg{0.0, 1.0};
  const std::vector<double> tau{0.0};
  const Vector2D idr(2, 1);
  const Vector2D lfcMismatched(1, 1);
  auto rpaInput = testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  EXPECT_THROW(thermoUtil::computeItcf(rpaInput, wvg, tau, 0.0, idr, lfcMismatched),
               std::runtime_error);
}

TEST(ThermoUtilTest, ComputeItcfInteractingPathSucceedsWithMatchingLfcShape) {
  auto in = testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  const std::vector<double> wvg{0.0, 1.0};
  const std::vector<double> tau{0.0, 0.3};
  Vector2D idr(2, 3);
  idr.fill(0, std::vector<double>{0.4, 0.2, 0.1});
  idr.fill(1, std::vector<double>{0.3, 0.15, 0.05});
  Vector2D lfc(2, 1);
  lfc(0, 0) = 0.1;
  lfc(1, 0) = 0.2;

  const Vector2D res = thermoUtil::computeItcf(in, wvg, tau, 0.0, idr, lfc);
  EXPECT_EQ(res.size(0), wvg.size());
  EXPECT_EQ(res.size(1), tau.size());
}
