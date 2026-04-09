#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "schemes/qstlsiet.hpp"

TEST(QstlsIetApiAndUtilTest, ConstructorRejectsGroundStateConfiguration) {
  auto ground = testFixtures::makeQstlsIetInput(
      "QSTLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  EXPECT_THROW((QstlsIet(ground)), std::runtime_error);
}

TEST(QstlsIetApiAndUtilTest, ConstructorAcceptsFiniteConfiguration) {
  auto finite = testFixtures::makeQstlsIetInput(
      "QSTLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.8, 2);
  EXPECT_NO_THROW((QstlsIet(finite)));
}

TEST(QstlsIetApiAndUtilTest, AdrIetReturnsZeroAtXZero) {
  const std::vector<double> wvg{0.0, 0.5, 1.0};
  const std::vector<double> ssf{0.0, 0.7, 1.0};
  const std::vector<double> lfc{0.0, 0.1, 0.2};
  const std::vector<double> bf{0.0, 0.0, 0.0};

  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  std::vector<std::shared_ptr<Interpolator1D>> lfci{
      std::make_shared<Interpolator1D>(wvg, lfc),
      std::make_shared<Interpolator1D>(wvg, lfc)};
  auto bfi = std::make_shared<Interpolator1D>(wvg, bf);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  Vector3D fixed(2, 3, 3);
  Vector2D res(3, 2);
  res.fill(0, 1.0);

  QstlsIetUtil::AdrIet adrIet(
      0.8, 0.0, 2.0, 0.0, ssfi, lfci, bfi, {}, itg2);
  adrIet.get(wvg, fixed, res);
  EXPECT_DOUBLE_EQ(res(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(res(0, 1), 0.0);
}

TEST(QstlsIetApiAndUtilTest, AdrFixedIetReturnsZeroAtXZero) {
  const std::vector<double> wvg{0.0, 0.5, 1.0};
  auto itg1 = std::make_shared<Integrator1D>(1.0e-8);
  Vector3D fixedIet(2, 3, 3);
  fixedIet.fill(1, 1.0);
  QstlsIetUtil::AdrFixedIet adrFixed(0.8, 0.0, 2.0, 0.0, 0.0, itg1);
  adrFixed.get(0, wvg, fixedIet);
  EXPECT_DOUBLE_EQ(fixedIet(0, 0, 0), 0.0);
}
