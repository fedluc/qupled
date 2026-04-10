#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "thermo/free_energy.hpp"

TEST(FreeEnergyTest, ComputesNonNormalizedIntegral) {
  const std::vector<double> rsGrid{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> rsu{2.0, 2.0, 2.0, 2.0};
  auto itp = std::make_shared<Interpolator1D>(rsGrid, rsu);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const FreeEnergy nonNormalized(2.0, itp, itg, false);
  EXPECT_NEAR(nonNormalized.get(), 4.0, 1.0e-6);
}

TEST(FreeEnergyTest, ComputesNormalizedIntegral) {
  const std::vector<double> rsGrid{0.0, 1.0, 2.0, 3.0};
  const std::vector<double> rsu{2.0, 2.0, 2.0, 2.0};
  auto itp = std::make_shared<Interpolator1D>(rsGrid, rsu);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const FreeEnergy normalized(2.0, itp, itg, true);
  EXPECT_NEAR(normalized.get(), 1.0, 1.0e-6);
}

TEST(FreeEnergyTest, ReturnsNegativeInfinityWhenNormalizedAtZeroCoupling) {
  const std::vector<double> rsGrid{0.0, 1.0, 2.0};
  const std::vector<double> rsu{1.0, 1.0, 1.0};
  auto itp = std::make_shared<Interpolator1D>(rsGrid, rsu);
  auto itg = std::make_shared<Integrator1D>(1.0e-8);

  const FreeEnergy normalizedAtZero(0.0, itp, itg, true);
  EXPECT_TRUE(std::isinf(normalizedAtZero.get()));
  EXPECT_LT(normalizedAtZero.get(), 0.0);
}
