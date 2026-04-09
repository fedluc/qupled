#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "thermo/itcf.hpp"

TEST(ItcfTest, NonInteractingFiniteReturnsExpectedXZeroFormulaIn2D) {
  auto in =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D2, 1.0, 0.7, 3);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  thermoUtil::ItcfNonInteracting itcf(
      0.0, in, 0.2, 0.0, 0.0, 2.0, itg, 0.0);

  const double expected = in->getDegeneracy() * (1.0 - std::exp(-1.0 / in->getDegeneracy()));
  EXPECT_NEAR(itcf.get(), expected, 1.0e-12);
}

TEST(ItcfTest, NonInteractingFiniteDelegatesTauZeroToStaticStructureFactor) {
  auto in =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  thermoUtil::ItcfNonInteracting itcf(
      0.0, in, 0.0, 0.0, 0.0, 2.0, itg, 0.0);

  const double v = itcf.get();
  EXPECT_TRUE(std::isfinite(v));
  EXPECT_GE(v, 0.0);
}

TEST(ItcfTest, NonInteractingGroundTauZeroMatchesSsfGroundReference) {
  auto in3d =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  thermoUtil::ItcfNonInteractingGround tauZero(1.0, in3d, 0.0, itg2);
  EXPECT_NEAR(tauZero.get(), 11.0 / 16.0, 1.0e-12);
}

TEST(ItcfTest, NonInteractingGroundReturnsZeroAtXZero) {
  auto in3d =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  thermoUtil::ItcfNonInteractingGround xZero(0.0, in3d, 0.2, itg2);
  EXPECT_DOUBLE_EQ(xZero.get(), 0.0);
}

TEST(ItcfTest, InteractingFiniteReturnsZeroAtXZero) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  const std::vector<double> idr{0.4, 0.2, 0.1};
  const std::vector<double> lfc{0.1};

  thermoUtil::Itcf xZero(0.0, in, 0.3, 0.5, lfc, idr);
  EXPECT_DOUBLE_EQ(xZero.get(), 0.0);
}

TEST(ItcfTest, InteractingFiniteReturnsHartreeFockValueAtZeroCoupling) {
  auto inZeroCoupling =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.7, 3);
  const std::vector<double> idr{0.4, 0.2, 0.1};
  const std::vector<double> lfc{0.1};
  thermoUtil::Itcf couplingZero(1.0, inZeroCoupling, 0.3, 0.5, lfc, idr);
  EXPECT_DOUBLE_EQ(couplingZero.get(), 0.5);
}

TEST(ItcfTest, InteractingFiniteTauZeroBranchReturnsFiniteValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  const std::vector<double> idr{0.4, 0.2, 0.1};
  const std::vector<double> lfc{0.1};
  thermoUtil::Itcf tauZero(1.0, in, 0.0, 0.8, lfc, idr);
  EXPECT_TRUE(std::isfinite(tauZero.get()));
}

TEST(ItcfTest, InteractingFiniteTauPositiveBranchReturnsFiniteValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  const std::vector<double> idr{0.4, 0.2, 0.1};
  const std::vector<double> lfc{0.1};
  thermoUtil::Itcf finiteTau(1.0, in, 0.3, 0.8, lfc, idr);
  EXPECT_TRUE(std::isfinite(finiteTau.get()));
}

TEST(ItcfTest, InteractingGroundReturnsZeroAtXZero) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  const std::vector<double> lfc{0.1};
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  thermoUtil::ItcfGround xZero(0.0, in, 0.3, 0.5, lfc, itg);
  EXPECT_DOUBLE_EQ(xZero.get(), 0.0);
}

TEST(ItcfTest, InteractingGroundReturnsHartreeFockValueAtZeroCoupling) {
  auto inZeroCoupling =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 0.0, 0.0, 2);
  const std::vector<double> lfc{0.1};
  auto itg = std::make_shared<Integrator1D>(1.0e-8);
  thermoUtil::ItcfGround couplingZero(1.0, inZeroCoupling, 0.3, 0.5, lfc, itg);
  EXPECT_DOUBLE_EQ(couplingZero.get(), 0.5);
}

TEST(ItcfTest, InteractingGroundTauZeroBranchReturnsFiniteValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  const std::vector<double> lfc{0.1};
  auto itg = std::make_shared<Integrator1D>(1.0e-8);
  thermoUtil::ItcfGround tauZero(1.0, in, 0.0, 0.8, lfc, itg);
  EXPECT_TRUE(std::isfinite(tauZero.get()));
}

TEST(ItcfTest, InteractingGroundTauPositiveBranchReturnsFiniteValue) {
  auto in =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  const std::vector<double> lfc{0.1};
  auto itg = std::make_shared<Integrator1D>(1.0e-8);
  thermoUtil::ItcfGround finiteTau(1.0, in, 0.3, 0.8, lfc, itg);
  EXPECT_TRUE(std::isfinite(finiteTau.get()));
}
