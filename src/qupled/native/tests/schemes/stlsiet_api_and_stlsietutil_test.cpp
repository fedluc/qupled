#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/stlsiet.hpp"

TEST(StlsIetApiAndUtilTest, ConstructorGuardsAndComputePath) {
  auto inBad = testFixtures::makeStlsIetInput(
      "STLS-IOI", "sqrt", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);
  EXPECT_THROW((StlsIet(inBad)), std::runtime_error);

  auto inOk = testFixtures::makeStlsIetInput(
      "STLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  StlsIet solver(inOk);
  EXPECT_EQ(solver.compute(), 0);
  EXPECT_EQ(solver.getBf().size(), solver.getWvg().size());
}

TEST(StlsIetApiAndUtilTest, SlfcUtilityCovers2D3DAndXZeroSpecialCase) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf{0.0, 0.4, 0.8, 1.0, 1.0};
  const std::vector<double> lfc{0.0, 0.1, 0.2, 0.2, 0.2};
  const std::vector<double> bf(wvg.size(), 0.0);
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto lfci = std::make_shared<Interpolator1D>(wvg, lfc);
  auto bfi = std::make_shared<Interpolator1D>(wvg, bf);
  auto itg = std::make_shared<Integrator2D>(1.0e-6);

  auto in3 = testFixtures::makeBaseInput(
      "STLS-HNC", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  auto in2 = testFixtures::makeBaseInput(
      "STLS-HNC", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);

  StlsIetUtil::Slfc x0_3d(0.0, 0.0, 2.0, ssfi, lfci, bfi, {}, itg, in3);
  StlsIetUtil::Slfc x1_3d(1.0, 0.0, 2.0, ssfi, lfci, bfi, {}, itg, in3);
  StlsIetUtil::Slfc x0_2d(0.0, 0.0, 2.0, ssfi, lfci, bfi, {}, itg, in2);

  EXPECT_DOUBLE_EQ(x0_3d.get(), 0.0);
  EXPECT_DOUBLE_EQ(x0_2d.get(), 0.0);
  EXPECT_TRUE(std::isfinite(x1_3d.get()));
}
