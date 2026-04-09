#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/rpa.hpp"

TEST(RpaApiAndUtilTest, ComputeFiniteStateBranch) {
  auto finite = testFixtures::makeBaseInput(
      "RPA", dimensionsUtil::Dimension::D3, 1.0, 0.8, 3);
  Rpa rpaFinite(finite, false);
  EXPECT_EQ(rpaFinite.compute(), 0);
  EXPECT_EQ(rpaFinite.getSsf().size(), rpaFinite.getWvg().size());
}

TEST(RpaApiAndUtilTest, ComputeGroundStateBranch) {
  auto ground =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  Rpa rpaGround(ground, false);
  EXPECT_EQ(rpaGround.compute(), 0);
  EXPECT_EQ(rpaGround.getSsf().size(), rpaGround.getWvg().size());
}

TEST(RpaApiAndUtilTest, UtilitySsfReturnsExpectedEdgeBehavior) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.7, 3);
  const std::vector<double> idr{0.3, 0.2, 0.1};
  const std::vector<double> lfcStatic{0.0};

  RpaUtil::Ssf s0(0.0, 0.9, lfcStatic, in, idr);
  EXPECT_DOUBLE_EQ(s0.get(), 0.0);
}

TEST(RpaApiAndUtilTest, UtilitySsfReturnsStaticBranchExpectedValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.7, 3);
  const std::vector<double> idr{0.3, 0.2, 0.1};
  const std::vector<double> lfcStatic{0.0};

  RpaUtil::Ssf sc(1.0, 0.9, lfcStatic, in, idr);
  EXPECT_NEAR(sc.get(), 0.9, 1.0e-14);
}

TEST(RpaApiAndUtilTest, UtilitySsfReturnsDynamicBranchExpectedValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.7, 3);
  const std::vector<double> idr{0.3, 0.2, 0.1};
  const std::vector<double> lfcDynamic{0.0, 0.0, 0.0};

  RpaUtil::Ssf sd(1.0, 0.9, lfcDynamic, in, idr);
  EXPECT_NEAR(sd.get(), 0.9, 1.0e-14);
}

TEST(RpaApiAndUtilTest, UtilitySsfGroundReturnsZeroAtXZero) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.0, 2);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);
  const std::vector<double> lfc{0.0};

  RpaUtil::SsfGround s0(0.0, 0.7, lfc, itg, in);
  EXPECT_DOUBLE_EQ(s0.get(), 0.0);
}

TEST(RpaApiAndUtilTest, UtilitySsfGroundUsesZeroCouplingFallback) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.0, 2);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);
  const std::vector<double> lfc{0.0};

  RpaUtil::SsfGround s1(1.0, 0.7, lfc, itg, in);
  EXPECT_NEAR(s1.get(), 0.7, 1.0e-14);
}
