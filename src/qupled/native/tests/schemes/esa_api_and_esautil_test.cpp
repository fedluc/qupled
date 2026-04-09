#include <gtest/gtest.h>

#include <cmath>

#include "fixtures/input_builders.hpp"
#include "schemes/esa.hpp"
#include "util/dual.hpp"

TEST(EsaApiAndUtilTest, EsaComputeProducesFiniteLocalFieldCorrection) {
  auto in =
      testFixtures::makeBaseInput("ESA", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  ESA esa(in);
  EXPECT_EQ(esa.compute(), 0);
  EXPECT_EQ(esa.getLfc().size(0), static_cast<int>(esa.getWvg().size()));
  EXPECT_TRUE(std::isfinite(esa.getLfc()(1, 0)));
}

TEST(EsaApiAndUtilTest, SlfcCachesCoefficientsAndExposesFiniteHelpers) {
  ESAUtil::Slfc slfc(1.0, 0.7);
  EXPECT_FALSE(slfc.coeff.valid);

  const double g1 = slfc.get(1.0);
  EXPECT_TRUE(slfc.coeff.valid);
  const double g2 = slfc.get(1.0);
  EXPECT_DOUBLE_EQ(g1, g2);

  EXPECT_TRUE(std::isfinite(slfc.nn(1.0)));
  EXPECT_TRUE(std::isfinite(slfc.csr(1.0)));
  EXPECT_TRUE(std::isfinite(slfc.activationFunction(1.0)));

  const auto fFinite = slfc.freeEnergy();
  EXPECT_TRUE(std::isfinite(fFinite.val()));

  ESAUtil::Slfc ground(1.0, 0.0);
  const auto fGround = ground.freeEnergy();
  EXPECT_TRUE(std::isfinite(fGround.val()));
}
