#include <gtest/gtest.h>

#include <cmath>

#include "util/num_util.hpp"

TEST(NumUtilTest, IsZeroUsesToleranceWindow) {
  EXPECT_TRUE(numUtil::isZero(0.0));
  EXPECT_TRUE(numUtil::isZero(0.5 * numUtil::dtol));
  EXPECT_FALSE(numUtil::isZero(2.0 * numUtil::dtol));
}

TEST(NumUtilTest, EqualTolTreatsCloseValuesAsEqual) {
  constexpr double x = 10.0;
  EXPECT_TRUE(numUtil::equalTol(x, x + 0.5 * x * numUtil::dtol));
  EXPECT_FALSE(numUtil::equalTol(x, x + 2.0 * x * numUtil::dtol));
}

TEST(NumUtilTest, LargerThanUsesRelativeTolerance) {
  constexpr double x = 10.0;
  constexpr double y_within_tol = x - 0.5 * x * numUtil::dtol;
  constexpr double y_beyond_tol = x - 2.0 * x * numUtil::dtol;

  EXPECT_FALSE(numUtil::largerThan(x, y_within_tol));
  EXPECT_TRUE(numUtil::largerThan(x, y_beyond_tol));
}

TEST(NumUtilTest, LambdaConstantsRemainSelfConsistent) {
  EXPECT_NEAR(numUtil::lambda * numUtil::lambda * numUtil::lambda,
              numUtil::lambda3,
              1.0e-15);
}
