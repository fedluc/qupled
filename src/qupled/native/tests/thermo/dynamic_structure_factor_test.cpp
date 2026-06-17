#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "thermo/dynamic_structure_factor.hpp"

namespace {

  constexpr double x = 1.2;
  constexpr double Omega = 0.7;
  constexpr double Theta = 0.8;
  constexpr double mu = 0.3;

  std::shared_ptr<Integrator1D> makeIntegrator() {
    return std::make_shared<Integrator1D>(Integrator1D::Type::SINGULAR, 1.0e-9);
  }

} // namespace

TEST(IdealDynamicResponseTest, DynamicRealPartMatchesReference) {
  thermoUtil::IdealDynamicResponse response(
      x, Omega, Theta, mu, makeIntegrator());
  EXPECT_NEAR(response.real(), 0.4245138149602567, 1.0e-8);
}

TEST(IdealDynamicResponseTest, DynamicImaginaryPartMatchesReference) {
  thermoUtil::IdealDynamicResponse response(
      x, Omega, Theta, mu, makeIntegrator());
  EXPECT_NEAR(response.imaginary(), 0.20032386853123246, 1.0e-12);
}

TEST(IdealDynamicResponseTest, StaticResponseMatchesReference) {
  thermoUtil::IdealDynamicResponse response(
      x, 0.0, Theta, mu, makeIntegrator());
  EXPECT_NEAR(response.real(), 0.4803089457348037, 1.0e-8);
  EXPECT_DOUBLE_EQ(response.imaginary(), 0.0);
}

TEST(IdealDynamicResponseTest, GridComputationUsesWaveVectorFrequencyLayout) {
  const std::vector<double> wvg{0.8, 1.2};
  const std::vector<double> frequency{0.0, 0.7};
  const auto [realPart, imaginaryPart] =
      thermoUtil::computeIdealDynamicResponse(
          wvg, frequency, Theta, mu, 1.0e-8);

  EXPECT_EQ(realPart.size(0), wvg.size());
  EXPECT_EQ(realPart.size(1), frequency.size());
  EXPECT_EQ(imaginaryPart.size(0), wvg.size());
  EXPECT_EQ(imaginaryPart.size(1), frequency.size());
  EXPECT_NEAR(realPart(1, 1), 0.4245138149602567, 1.0e-7);
  EXPECT_NEAR(imaginaryPart(1, 1), 0.20032386853123246, 1.0e-12);
}

TEST(IdealDynamicResponseTest, RejectsUnsupportedDomain) {
  EXPECT_THROW(
      thermoUtil::IdealDynamicResponse(0.0, Omega, Theta, mu, makeIntegrator()),
      std::runtime_error);
  EXPECT_THROW(
      thermoUtil::IdealDynamicResponse(x, -Omega, Theta, mu, makeIntegrator()),
      std::runtime_error);
  EXPECT_THROW(
      thermoUtil::IdealDynamicResponse(x, Omega, 0.0, mu, makeIntegrator()),
      std::runtime_error);
  EXPECT_THROW(
      thermoUtil::computeIdealDynamicResponse({x}, {Omega}, Theta, mu, 0.0),
      std::runtime_error);
}

TEST(DynamicStructureFactorTest, MatchesDynamicAndStaticReferences) {
  auto in = testFixtures::makeBaseInput(
      "STLS", dimensionsUtil::Dimension::D3, 2.0, Theta, 4);
  const std::vector<double> wvg{x};
  const std::vector<double> frequency{0.0, Omega};
  const Vector2D lfc(std::vector<std::vector<double>>{{0.25}});

  const Vector2D result =
      thermoUtil::computeDSF(in, wvg, frequency, mu, lfc);

  EXPECT_NEAR(result(0, 0), 0.04345726902387644, 1.0e-9);
  EXPECT_NEAR(result(0, 1), 0.06462738121282116, 1.0e-9);
}

TEST(DynamicStructureFactorTest, RejectsUnsupportedInputs) {
  auto in2D = testFixtures::makeBaseInput(
      "STLS", dimensionsUtil::Dimension::D2, 2.0, Theta, 4);
  auto inGround = testFixtures::makeBaseInput(
      "STLS", dimensionsUtil::Dimension::D3, 2.0, 0.0, 4);
  auto inFinite = testFixtures::makeBaseInput(
      "STLS", dimensionsUtil::Dimension::D3, 2.0, Theta, 4);
  const Vector2D staticLfc(std::vector<std::vector<double>>{{0.25}});
  const Vector2D dynamicLfc(std::vector<std::vector<double>>{{0.25, 0.2}});

  EXPECT_THROW(thermoUtil::computeDSF(
                   in2D, {x}, {Omega}, mu, staticLfc),
               std::runtime_error);
  EXPECT_THROW(thermoUtil::computeDSF(
                   inGround, {x}, {Omega}, mu, staticLfc),
               std::runtime_error);
  EXPECT_THROW(thermoUtil::computeDSF(
                   inFinite, {x}, {Omega}, mu, dynamicLfc),
               std::runtime_error);
  EXPECT_THROW(thermoUtil::computeDSF(
                   inFinite, {0.0}, {-Omega}, mu, staticLfc),
               std::runtime_error);
}

TEST(DynamicStructureFactorTest, ReturnsZeroAtZeroWaveVector) {
  auto in = testFixtures::makeBaseInput(
      "STLS", dimensionsUtil::Dimension::D3, 2.0, Theta, 4);
  const Vector2D lfc(std::vector<std::vector<double>>{{0.25}});

  const Vector2D result = thermoUtil::computeDSF(
      in, {0.0}, {0.0, Omega}, mu, lfc);

  EXPECT_DOUBLE_EQ(result(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(result(0, 1), 0.0);
}

TEST(DynamicStructureFactorTransformTest,
     AdaptiveTransformUsesCquadInterpolatedBaseline) {
  auto in = testFixtures::makeBaseInput(
      "RPA", dimensionsUtil::Dimension::D3, 2.0, 1.0, 4);
  const std::vector<double> frequency{0.0, 1.0, 2.0};
  const Vector2D dsf(std::vector<std::vector<double>>{{0.0, 1.0, 2.0}});
  const Vector2D lfc(std::vector<std::vector<double>>{{0.0}});
  const double expected = 12.0 * (1.0 - 2.0 / std::exp(1.0));

  const Vector2D result =
      thermoUtil::computeAdaptiveITCFfromDSF(
          in, {4.0}, frequency, {0.5}, mu, lfc, dsf);

  EXPECT_NEAR(result(0, 0), expected, 1.0e-8);
}
