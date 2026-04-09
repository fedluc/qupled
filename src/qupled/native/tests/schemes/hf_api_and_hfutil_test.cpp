#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/hf.hpp"

TEST(HfApiAndUtilTest, ConstructorInitializesGridAndStorageShapes) {
  auto in =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.8, 3);

  HF hf(in, false);

  EXPECT_FALSE(hf.getWvg().empty());
  EXPECT_EQ(hf.getSsf().size(), hf.getWvg().size());
  EXPECT_EQ(hf.getIdr().size(0), hf.getWvg().size());
  EXPECT_EQ(hf.getIdr().size(1), static_cast<size_t>(in->getNMatsubara()));
  EXPECT_EQ(hf.getLfc().size(0), hf.getWvg().size());
  EXPECT_EQ(hf.getLfc().size(1), 1u);
}

TEST(HfApiAndUtilTest, ComputeFiniteStateExposesExpectedShapes) {
  auto finite = testFixtures::makeBaseInput(
      "HF", dimensionsUtil::Dimension::D3, 1.0, 0.8, 3);
  HF hfFinite(finite, false);
  EXPECT_EQ(hfFinite.compute(), 0);
  EXPECT_FALSE(hfFinite.getWvg().empty());
  EXPECT_EQ(hfFinite.getSsf().size(), hfFinite.getWvg().size());
  EXPECT_EQ(hfFinite.getIdr().size(0), static_cast<int>(hfFinite.getWvg().size()));
  EXPECT_FALSE(hfFinite.getSdr().empty());
  EXPECT_TRUE(std::isfinite(hfFinite.getUInt()));
  EXPECT_TRUE(std::isfinite(hfFinite.getChemicalPotential()));
}

TEST(HfApiAndUtilTest, ComputeGroundStateExposesExpectedShapes) {
  auto ground =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D2, 1.0, 0.0, 2);
  HF hfGround(ground, false);
  EXPECT_EQ(hfGround.compute(), 0);
  EXPECT_EQ(hfGround.getSsf().size(), hfGround.getWvg().size());
  EXPECT_TRUE(hfGround.getSdr().empty());
}

TEST(HfApiAndUtilTest, IdrGroundMatchesKnownAnalyticLimits) {
  EXPECT_DOUBLE_EQ(HFUtil::IdrGround(dimensionsUtil::Dimension::D3, 0.0, 0.0).get(), 1.0);
  EXPECT_DOUBLE_EQ(HFUtil::IdrGround(dimensionsUtil::Dimension::D3, 0.0, 1.0).get(), 0.0);
  EXPECT_DOUBLE_EQ(HFUtil::IdrGround(dimensionsUtil::Dimension::D2, 0.0, 0.0).get(),
                   2.0 / 3.0);
}

TEST(HfApiAndUtilTest, SsfGroundMatchesKnownAnalyticLimits) {
  EXPECT_NEAR(HFUtil::SsfGround(dimensionsUtil::Dimension::D3, 1.0).get(), 11.0 / 16.0,
              1.0e-12);
  EXPECT_DOUBLE_EQ(HFUtil::SsfGround(dimensionsUtil::Dimension::D3, 3.0).get(), 1.0);
  EXPECT_DOUBLE_EQ(HFUtil::SsfGround(dimensionsUtil::Dimension::D2, 0.0).get(), 0.0);
  EXPECT_DOUBLE_EQ(HFUtil::SsfGround(dimensionsUtil::Dimension::D2, 3.0).get(), 1.0);
}

TEST(HfApiAndUtilTest, FiniteTemperatureUtilityObjectsReturnFiniteValues) {
  auto in =
      testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  auto itg1 = std::make_shared<Integrator1D>(1.0e-8);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  HFUtil::Idr idr(in, 0.5, 0.0, 0.0, 2.0, itg1);
  const auto idrVals = idr.get();
  ASSERT_EQ(idrVals.size(), static_cast<size_t>(in->getNMatsubara()));
  for (const double v : idrVals) {
    EXPECT_TRUE(std::isfinite(v));
  }

  HFUtil::Ssf ssf(in, 0.5, 0.0, 0.0, 2.0, itg1, {}, itg2, idrVals[0]);
  EXPECT_TRUE(std::isfinite(ssf.get()));
}
