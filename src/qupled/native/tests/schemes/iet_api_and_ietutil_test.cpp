#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/iet.hpp"

TEST(IetApiAndUtilTest, BridgeFunctionReturnsZeroForHncTheory) {
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  IetUtil::BridgeFunction hnc("STLS-HNC", "sqrt", 1.0, 0.7, 1.0, itg);
  EXPECT_DOUBLE_EQ(hnc.get(), 0.0);
}

TEST(IetApiAndUtilTest, BridgeFunctionRejectsUnknownTheory) {
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  IetUtil::BridgeFunction invalidTheory("BAD", "sqrt", 1.0, 0.7, 1.0, itg);
  EXPECT_THROW(invalidTheory.get(), std::runtime_error);
}

TEST(IetApiAndUtilTest, BridgeFunctionRejectsStandardMappingAtGroundState) {
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  IetUtil::BridgeFunction badStandardGround(
      "STLS-IOI", "standard", 1.0, 0.0, 1.0, itg);
  EXPECT_THROW(badStandardGround.get(), std::runtime_error);
}

TEST(IetApiAndUtilTest, BridgeFunctionRejectsIoiOutsideSupportedRange) {
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  IetUtil::BridgeFunction ioiOutOfRange("STLS-IOI", "sqrt", 0.2, 1.0, 1.0, itg);
  EXPECT_THROW(ioiOutOfRange.get(), std::runtime_error);
}

TEST(IetApiAndUtilTest, BridgeFunctionRejectsLctOutsideSupportedRange) {
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  IetUtil::BridgeFunction lctOutOfRange("STLS-LCT", "sqrt", 0.2, 1.0, 1.0, itg);
  EXPECT_THROW(lctOutOfRange.get(), std::runtime_error);
}

TEST(IetApiAndUtilTest, IetInitProducesBridgeFunctionGrid) {
  auto in = testFixtures::makeStlsIetInput(
      "STLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Guess guess;
  guess.wvg = {0.0, 0.5, 1.0};
  guess.ssf = {0.0, 0.6, 1.0};
  guess.lfc = Vector2D(3, 2);
  for (int i = 0; i < 3; ++i) {
    guess.lfc(i, 0) = 0.1 * i;
    guess.lfc(i, 1) = 0.2 * i;
  }
  in->setGuess(guess);

  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5};
  Iet iet(in, wvg);
  iet.init();
  EXPECT_EQ(iet.getBf().size(), wvg.size());
}

TEST(IetApiAndUtilTest, InitialGuessFromInputInterpolatesToSolverGrid) {
  auto in = testFixtures::makeStlsIetInput(
      "STLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Guess guess;
  guess.wvg = {0.0, 0.5, 1.0};
  guess.ssf = {0.0, 0.6, 1.0};
  guess.lfc = Vector2D(3, 2);
  for (int i = 0; i < 3; ++i) {
    guess.lfc(i, 0) = 0.1 * i;
    guess.lfc(i, 1) = 0.2 * i;
  }
  in->setGuess(guess);

  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5};
  Iet iet(in, wvg);

  Vector2D lfc(static_cast<int>(wvg.size()), 2);
  EXPECT_TRUE(iet.initialGuessFromInput(lfc));
  EXPECT_DOUBLE_EQ(lfc(3, 0), 0.0);
  EXPECT_DOUBLE_EQ(lfc(3, 1), 0.0);
}
