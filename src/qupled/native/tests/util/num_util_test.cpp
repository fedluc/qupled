#include <gtest/gtest.h>

#include <cmath>

#include "util/num_util.hpp"

TEST(NumUtilTest, IsZeroReturnsTrueAtExactlyZero) {
  EXPECT_TRUE(numUtil::isZero(0.0));
}

TEST(NumUtilTest, IsZeroReturnsTrueForPositiveValueInsideTolerance) {
  EXPECT_TRUE(numUtil::isZero(0.5 * numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsTrueForNegativeValueInsideTolerance) {
  EXPECT_TRUE(numUtil::isZero(-0.5 * numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsFalseAtPositiveToleranceBoundary) {
  EXPECT_FALSE(numUtil::isZero(numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsFalseAtNegativeToleranceBoundary) {
  EXPECT_FALSE(numUtil::isZero(-numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsFalseForPositiveValueOutsideTolerance) {
  EXPECT_FALSE(numUtil::isZero(2.0 * numUtil::dtol));
}

TEST(NumUtilTest, IsZeroReturnsFalseForNegativeValueOutsideTolerance) {
  EXPECT_FALSE(numUtil::isZero(-2.0 * numUtil::dtol));
}

TEST(EqualTolTest, ReturnsTrueWhenDifferenceIsBelowRelativeTolerance) {
  constexpr double x = 10.0;
  EXPECT_TRUE(numUtil::equalTol(x, x + 0.5 * x * numUtil::dtol));
}

TEST(EqualTolTest, ReturnsFalseAtRelativeToleranceBoundary) {
  constexpr double x = 10.0;
  EXPECT_FALSE(numUtil::equalTol(x, x + x * numUtil::dtol));
}

TEST(EqualTolTest, ReturnsFalseWhenDifferenceExceedsRelativeTolerance) {
  constexpr double x = 10.0;
  EXPECT_FALSE(numUtil::equalTol(x, x + 2.0 * x * numUtil::dtol));
}

TEST(LargerThanTest, ReturnsFalseWhenGapIsBelowRelativeTolerance) {
  constexpr double x = 10.0;
  EXPECT_FALSE(numUtil::largerThan(x, x - 0.5 * x * numUtil::dtol));
}

TEST(LargerThanTest, ReturnsTrueWhenGapExceedsRelativeTolerance) {
  constexpr double x = 10.0;
  EXPECT_TRUE(numUtil::largerThan(x, x - 2.0 * x * numUtil::dtol));
}

TEST(NumUtilConstantsTest, InfIsInfinite) {
  EXPECT_TRUE(std::isinf(numUtil::Inf));
}

TEST(NumUtilConstantsTest, NaNIsNan) { EXPECT_TRUE(std::isnan(numUtil::NaN)); }

TEST(NumUtilConstantsTest, LambdaCubedMatchesLambda3Constant) {
  EXPECT_NEAR(numUtil::lambda * numUtil::lambda * numUtil::lambda,
              numUtil::lambda3,
              1.0e-15);
}
