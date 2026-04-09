#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "fixtures/tolerance.hpp"
#include "thermo/thermo_util.hpp"

TEST(ThermoUtilTest, FreeEnergyThrowsIfCouplingExceedsGridRange) {
  const std::vector<double> grid{0.0, 1.0, 2.0};
  const std::vector<double> rsu{0.0, 1.0, 2.0};

  EXPECT_THROW(thermoUtil::computeFreeEnergy(grid, rsu, 2.1), std::runtime_error);
}

TEST(ThermoUtilTest, ComputeItcfReturnsEmptyWhenInputsAreEmpty) {
  auto in = testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  const Vector2D idr(1, 1);
  const Vector2D lfc(1, 1);

  const Vector2D empty_wvg =
      thermoUtil::computeItcf(in, {}, {0.0}, 0.0, idr, lfc);
  EXPECT_TRUE(empty_wvg.empty());

  const Vector2D empty_tau =
      thermoUtil::computeItcf(in, {0.0, 1.0}, {}, 0.0, idr, lfc);
  EXPECT_TRUE(empty_tau.empty());
}

TEST(ThermoUtilTest, ComputeItcfReturnsHartreeFockBranchWithoutLfcValidation) {
  auto in = testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  const std::vector<double> wvg{0.0, 1.0};
  const std::vector<double> tau{0.0};
  const Vector2D idr(2, 1);
  const Vector2D lfc_mismatched(1, 1);

  const Vector2D result = thermoUtil::computeItcf(in, wvg, tau, 0.0, idr, lfc_mismatched);
  ASSERT_EQ(result.size(0), 2u);
  ASSERT_EQ(result.size(1), 1u);
  EXPECT_NEAR(result(0, 0), 0.0, testFixtures::kTightTol);
  EXPECT_NEAR(result(1, 0), 11.0 / 16.0, testFixtures::kTightTol);
}

TEST(ThermoUtilTest, ComputeItcfValidatesLfcRowsForInteractingTheories) {
  auto in = testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 4);
  const std::vector<double> wvg{0.0, 1.0};
  const std::vector<double> tau{0.0};
  const Vector2D idr(2, 1);
  const Vector2D lfc_mismatched(1, 1);

  EXPECT_THROW(thermoUtil::computeItcf(in, wvg, tau, 0.0, idr, lfc_mismatched),
               std::runtime_error);
}
