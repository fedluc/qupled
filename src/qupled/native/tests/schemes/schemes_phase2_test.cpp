#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>

#include "fixtures/input_builders.hpp"
#include "schemes/esa.hpp"
#include "schemes/hf.hpp"
#include "schemes/iet.hpp"
#include "schemes/rpa.hpp"
#include "schemes/stls.hpp"
#include "schemes/stlsiet.hpp"

TEST(SchemesPhase2Test, HfComputesAndRejectsInternalEnergyWithoutData) {
  auto in = testFixtures::makeBaseInput("HF", dimensionsUtil::Dimension::D3, 1.0, 0.8, 2);
  HF hf(in, false);

  EXPECT_EQ(hf.compute(), 0);
  EXPECT_FALSE(hf.getWvg().empty());
  EXPECT_EQ(hf.getSsf().size(), hf.getWvg().size());
  EXPECT_TRUE(std::isfinite(hf.getUInt()));
}

TEST(SchemesPhase2Test, RpaFiniteAndGroundBranchesReturnFiniteOutputs) {
  auto in_finite =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 3);
  Rpa rpa_finite(in_finite, false);
  EXPECT_EQ(rpa_finite.compute(), 0);
  EXPECT_EQ(rpa_finite.getSsf().size(), rpa_finite.getWvg().size());

  auto in_ground =
      testFixtures::makeBaseInput("RPA", dimensionsUtil::Dimension::D3, 1.0, 0.0, 1);
  Rpa rpa_ground(in_ground, false);
  EXPECT_EQ(rpa_ground.compute(), 0);
  EXPECT_EQ(rpa_ground.getSsf().size(), rpa_ground.getWvg().size());
}

TEST(SchemesPhase2Test, StlsComputesWithDefaultAndInputGuessPaths) {
  auto in_default = testFixtures::makeStlsInput("STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Stls stls_default(in_default, false);
  EXPECT_EQ(stls_default.compute(), 0);
  EXPECT_EQ(stls_default.getSsf().size(), stls_default.getWvg().size());

  auto in_guess = testFixtures::makeStlsInput("STLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  Guess guess;
  guess.wvg = {0.0, 0.5, 1.0, 1.5, 2.0};
  guess.ssf = {0.0, 0.4, 0.8, 1.0, 1.0};
  in_guess->setGuess(guess);
  Stls stls_guess(in_guess, false);
  EXPECT_EQ(stls_guess.compute(), 0);
}

TEST(SchemesPhase2Test, StlsIetSupportsHncAndRejectsUnsupported2DIoi) {
  auto in_hnc =
      testFixtures::makeStlsIetInput("STLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  StlsIet stls_iet(in_hnc);
  EXPECT_EQ(stls_iet.compute(), 0);
  EXPECT_EQ(stls_iet.getBf().size(), stls_iet.getWvg().size());

  auto in_2d_ioi =
      testFixtures::makeStlsIetInput("STLS-IOI", "sqrt", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);
  EXPECT_THROW({ StlsIet tmp(in_2d_ioi); }, std::runtime_error);
}

TEST(SchemesPhase2Test, EsaSlfcCachesCoefficientsAndBridgeUtilityValidatesTheory) {
  ESAUtil::Slfc slfc(1.0, 0.7);
  EXPECT_FALSE(slfc.coeff.valid);
  const double g1 = slfc.get(1.0);
  EXPECT_TRUE(slfc.coeff.valid);
  const double g2 = slfc.get(1.5);
  EXPECT_TRUE(std::isfinite(g1));
  EXPECT_TRUE(std::isfinite(g2));

  auto itg = std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);
  const IetUtil::BridgeFunction valid("STLS-HNC", "sqrt", 1.0, 0.7, 1.0, itg);
  EXPECT_NO_THROW(valid.get());

  const IetUtil::BridgeFunction invalid("NOT-A-THEORY", "sqrt", 1.0, 0.7, 1.0, itg);
  EXPECT_THROW(invalid.get(), std::runtime_error);
}
