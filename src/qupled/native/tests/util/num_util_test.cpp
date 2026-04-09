#include <gtest/gtest.h>

#include <cmath>

#include "util/num_util.hpp"

TEST(NumUtilTest, IsZeroReturnsTrueForValuesWithinTolerance) {
  EXPECT_TRUE(numUtil::isZero(0.0));
  EXPECT_TRUE(numUtil::isZero(0.5 * numUtil::dtol));
  EXPECT_TRUE(numUtil::isZero(-0.5 * numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsFalseForValuesOutsideTolerance) {
  EXPECT_FALSE(numUtil::isZero(2.0 * numUtil::dtol));
  EXPECT_FALSE(numUtil::isZero(-2.0 * numUtil::dtol));
}

TEST(NumUtilTest, EqualTolReturnsTrueWhenDifferenceIsBelowThreshold) {
  constexpr double x = 10.0;
  EXPECT_TRUE(numUtil::equalTol(x, x + 0.5 * x * numUtil::dtol));
}

TEST(NumUtilTest, EqualTolReturnsFalseAtOrAboveThreshold) {
  constexpr double x = 10.0;
  EXPECT_FALSE(numUtil::equalTol(x, x + x * numUtil::dtol));
  EXPECT_FALSE(numUtil::equalTol(x, x + 2.0 * x * numUtil::dtol));
}

TEST(NumUtilTest, LargerThanReturnsTrueOnlyBeyondRelativeTolerance) {
  constexpr double x = 10.0;
  EXPECT_FALSE(numUtil::largerThan(x, x - 0.5 * x * numUtil::dtol));
  EXPECT_TRUE(numUtil::largerThan(x, x - 2.0 * x * numUtil::dtol));
}

TEST(NumUtilTest, ExposesExpectedSentinelAndLambdaConstants) {
  EXPECT_TRUE(std::isinf(numUtil::Inf));
  EXPECT_TRUE(std::isnan(numUtil::NaN));
  EXPECT_NEAR(numUtil::lambda * numUtil::lambda * numUtil::lambda,
              numUtil::lambda3,
              1.0e-15);
}
