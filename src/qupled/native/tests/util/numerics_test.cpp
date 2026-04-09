#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "fixtures/tolerance.hpp"
#include "util/numerics.hpp"

TEST(NumericsTest, SpecialFunctionsHandleBasicReferenceValues) {
  EXPECT_TRUE(std::isinf(SpecialFunctions::coth(0.0)));
  EXPECT_NEAR(SpecialFunctions::besselJ0(0.0), 1.0, 1.0e-12);
  EXPECT_TRUE(std::isinf(SpecialFunctions::ellipticK(1.0)));
  EXPECT_TRUE(std::isinf(SpecialFunctions::ellipticE(1.0)));
  EXPECT_GT(SpecialFunctions::fermiDirac12(1.0), SpecialFunctions::fermiDirac12(0.0));
  EXPECT_GT(SpecialFunctions::fermiDiracm12(1.0),
            SpecialFunctions::fermiDiracm12(0.0));
}

TEST(NumericsTest, Interpolator1DSupportsEvaluationClampingAndReset) {
  const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> y{0.0, 1.0, 4.0, 9.0};

  Interpolator1D itp(x, y);
  EXPECT_TRUE(itp.isValid());
  EXPECT_NEAR(itp.eval(2.0), 4.0, 1.0e-12);
  EXPECT_NEAR(itp.eval(5.0), 9.0, 1.0e-12);

  const std::vector<double> y_reset{0.0, 1.0, 8.0, 27.0};
  itp.reset(x[0], y_reset[0], x.size());
  EXPECT_TRUE(itp.isValid());
  EXPECT_NEAR(itp.eval(2.0), 8.0, 1.0e-12);
}

TEST(NumericsTest, Interpolator2DInterpolatesRegularGridValues) {
  const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> y{0.0, 1.0, 2.0, 3.0};
  std::vector<double> z;
  z.reserve(x.size() * y.size());
  for (double xi : x) {
    for (double yi : y) {
      z.push_back(xi + yi);
    }
  }

  Interpolator2D itp(
      x[0], y[0], z[0], static_cast<int>(x.size()), static_cast<int>(y.size()));
  EXPECT_TRUE(itp.isValid());
  EXPECT_NEAR(itp.eval(1.5, 1.5), 3.0, 2.0e-2);
  EXPECT_NEAR(itp.eval(3.0, 0.0), 3.0, 2.0e-2);

  Interpolator2D itp_invalid(x[0], y[0], z[0], 1, 1);
  EXPECT_FALSE(itp_invalid.isValid());
}

TEST(NumericsTest, RootSolversConvergeForSimpleNonlinearEquation) {
  const auto func = [](double x) { return x * x - 2.0; };

  BrentRootSolver brent;
  brent.solve(func, {1.0, 2.0});
  EXPECT_NEAR(brent.getSolution(), std::sqrt(2.0), 1.0e-9);

  SecantSolver secant;
  secant.solve(func, {1.0, 2.0});
  EXPECT_NEAR(secant.getSolution(), std::sqrt(2.0), 1.0e-9);
}

TEST(NumericsTest, RootSolversReportInvalidOrNonConvergedSetups) {
  BrentRootSolver brent;
  EXPECT_THROW(brent.solve([](double x) { return x * x + 1.0; }, {-1.0, 1.0}),
               std::runtime_error);

  SecantSolver secant(1.0e-12, 2);
  EXPECT_THROW(secant.solve([](double x) { return x * x + 1.0; }, {0.0, 1.0}),
               std::runtime_error);
}

TEST(NumericsTest, Integrator1DComputesSimpleReferenceIntegrals) {
  Integrator1D itg_default(1.0e-8);
  itg_default.compute([](double x) { return x; }, Integrator1D::Param(0.0, 1.0));
  EXPECT_NEAR(itg_default.getSolution(), 0.5, 1.0e-6);

  Integrator1D itg_singular(Integrator1D::Type::SINGULAR, 1.0e-8);
  itg_singular.compute(
      [](double x) { return 1.0 / std::sqrt(x); }, Integrator1D::Param(0.0, 1.0));
  EXPECT_NEAR(itg_singular.getSolution(), 2.0, 1.0e-5);

  Integrator1D itg_fourier(Integrator1D::Type::FOURIER, 1.0e-6);
  EXPECT_THROW(itg_fourier.compute([](double x) { return x; }, Integrator1D::Param(0.0, 1.0)),
               std::runtime_error);
}

TEST(NumericsTest, Integrator2DComputesSeparableIntegralWithAndWithoutGrid) {
  Integrator2D itg(1.0e-6);
  const auto fx = [](double x) { return x; };
  const auto fy = [](double y) { return y; };
  const Integrator2D::Param param(0.0, 1.0, 0.0, 1.0);

  itg.compute(fx, fy, param, {});
  EXPECT_NEAR(itg.getSolution(), 0.25, testFixtures::kLooseTol);

  itg.compute(fx, fy, param, {0.0, 0.25, 0.5, 0.75, 1.0});
  EXPECT_NEAR(itg.getSolution(), 0.25, 2.0e-3);
}
