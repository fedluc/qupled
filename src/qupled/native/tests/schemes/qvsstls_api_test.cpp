#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/qvsstls.hpp"
#include "vs/grid_point.hpp"

TEST(QvsStlsApiTest, ConstructorGuardsAndManagerContracts) {
  auto ground = testFixtures::makeQVSStlsInput(
      "QVSSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  EXPECT_THROW((QVSStls(ground)), std::runtime_error);

  auto finite = testFixtures::makeQVSStlsInput(
      "QVSSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.8, 2);
  VSQstlsManager mgr(finite);
  mgr.setAlpha(0.4);
  EXPECT_DOUBLE_EQ(mgr.getAlpha(), 0.4);
  EXPECT_FALSE(static_cast<const VSManager &>(mgr).getWvg(GridPoints::CENTER).empty());
}

TEST(QvsStlsApiTest, QAdderKeepsUnitStructureFactorBaselineNearZero) {
  const std::vector<double> wvg{0.1, 0.4, 0.8, 1.2, 1.6};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg1 = std::make_shared<Integrator1D>(1.0e-8);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  const QAdder qadder(0.7, 0.0, wvg.front(), wvg.back(), {}, itg1, itg2, ssfi);
  EXPECT_NEAR(qadder.get(), 0.0, 1.0e-6);
}
