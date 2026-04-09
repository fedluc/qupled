#include <gtest/gtest.h>

#include "util/dual.hpp"

TEST(DualTest, SeedsOnlyTheRequestedGradientEntry) {
  Dual<1> x(2.0, 3, 1);

  ASSERT_EQ(x.grad.size(), 3u);
  EXPECT_DOUBLE_EQ(x.func, 2.0);
  EXPECT_DOUBLE_EQ(x.grad[0], 0.0);
  EXPECT_DOUBLE_EQ(x.grad[1], 1.0);
  EXPECT_DOUBLE_EQ(x.grad[2], 0.0);
}

TEST(DualTest, MultiplicationPropagatesFirstDerivative) {
  Dual<1> x(3.0, 1, 0);
  const auto y = x * x;

  EXPECT_DOUBLE_EQ(y.func, 9.0);
  ASSERT_EQ(y.grad.size(), 1u);
  EXPECT_DOUBLE_EQ(y.grad[0], 6.0);
}
