#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/stls.hpp"

TEST(StlsApiAndUtilTest, ConstructorAcceptsValidInput) {
  auto in =
      testFixtures::makeStlsInput("STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);

  EXPECT_NO_THROW((Stls(in, false)));
}

TEST(StlsApiAndUtilTest, ComputeSupportsDefaultAndUserProvidedGuess) {
  auto inDefault = testFixtures::makeStlsInput(
      "STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Stls stlsDefault(inDefault, false);
  EXPECT_EQ(stlsDefault.compute(), 0);
  EXPECT_EQ(stlsDefault.getSsf().size(), stlsDefault.getWvg().size());
  EXPECT_GE(stlsDefault.getError(), 0.0);

  auto inGuess = testFixtures::makeStlsInput(
      "STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Guess guess;
  guess.wvg = {0.0, 0.5, 1.0, 1.5, 2.0};
  guess.ssf = {0.0, 0.4, 0.8, 1.0, 1.0};
  inGuess->setGuess(guess);
  Stls stlsGuess(inGuess, false);
  EXPECT_EQ(stlsGuess.compute(), 0);
}

TEST(StlsApiAndUtilTest, DynamicPointerCastHelpersCoverSuccessAndFailure) {
  auto stlsIn = std::make_shared<StlsInput>();
  stlsIn->setTheory("STLS");
  std::shared_ptr<const Input> asBase = stlsIn;
  auto castOk = StlsUtil::dynamic_pointer_cast<Input, StlsInput>(asBase);
  EXPECT_EQ(castOk.get(), stlsIn.get());

  auto plainInput = std::make_shared<Input>();
  plainInput->setTheory("HF");
  std::shared_ptr<const Input> plainBase = plainInput;
  EXPECT_THROW((StlsUtil::dynamic_pointer_cast<Input, StlsInput>(plainBase)),
               std::runtime_error);
}

TEST(StlsApiAndUtilTest, SlfcUtilityHandlesSpecialPointsIn2DAnd3D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  auto in3 =
      testFixtures::makeBaseInput("STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  auto in2 =
      testFixtures::makeBaseInput("STLS", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);

  StlsUtil::Slfc g0_3d(0.0, 0.0, 2.0, ssfi, itg, in3);
  StlsUtil::Slfc gx_3d(1.0, 0.0, 2.0, ssfi, itg, in3);
  StlsUtil::Slfc g0_2d(0.0, 0.0, 2.0, ssfi, itg, in2);
  EXPECT_DOUBLE_EQ(g0_3d.get(), 0.0);
  EXPECT_DOUBLE_EQ(g0_2d.get(), 0.0);
  EXPECT_TRUE(std::isfinite(gx_3d.get()));
}
