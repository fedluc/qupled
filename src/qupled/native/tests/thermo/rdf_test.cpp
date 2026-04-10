#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "thermo/rdf.hpp"

TEST(RdfTest, StaysAtOneForUnitStructureFactorIn2D) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-8);
  auto itgf =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  Rdf rdf2d(0.3, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D2);
  EXPECT_NEAR(rdf2d.get(), 1.0, 1.0e-6);
}

TEST(RdfTest, StaysAtOneForUnitStructureFactorIn3DAtZeroDistance) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-8);
  auto itgf =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  Rdf rdf3dAtZero(
      0.0, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(rdf3dAtZero.get(), 1.0, 1.0e-6);
}

TEST(RdfTest, StaysAtOneForUnitStructureFactorIn3DAtFiniteDistance) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto itp = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg =
      std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT, 1.0e-8);
  auto itgf =
      std::make_shared<Integrator1D>(Integrator1D::Type::FOURIER, 1.0e-8);

  Rdf rdf3d(1.2, wvg.back(), itp, itg, itgf, dimensionsUtil::Dimension::D3);
  EXPECT_NEAR(rdf3d.get(), 1.0, 1.0e-5);
}
