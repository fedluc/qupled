#include <gtest/gtest.h>

#include <cmath>

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

TEST(DualTest, ArithmeticWithScalarsKeepsExpectedGradient) {
  const Dual<1> x(4.0, 1, 0);
  const auto y = (2.0 * x + 3.0) / 5.0;

  EXPECT_DOUBLE_EQ(y.func, 2.2);
  ASSERT_EQ(y.grad.size(), 1u);
  EXPECT_DOUBLE_EQ(y.grad[0], 0.4);
}

TEST(DualTest, ElementaryFunctionsFollowChainRule) {
  const Dual<1> x(2.0, 1, 0);

  const auto e = exp(x);
  EXPECT_DOUBLE_EQ(e.func, std::exp(2.0));
  EXPECT_DOUBLE_EQ(e.grad[0], std::exp(2.0));

  const auto l = log(x);
  EXPECT_DOUBLE_EQ(l.func, std::log(2.0));
  EXPECT_DOUBLE_EQ(l.grad[0], 0.5);

  const auto s = sqrt(x);
  EXPECT_DOUBLE_EQ(s.func, std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(s.grad[0], 0.5 / std::sqrt(2.0));

  const auto t = tanh(x);
  EXPECT_DOUBLE_EQ(t.func, std::tanh(2.0));
  EXPECT_DOUBLE_EQ(t.grad[0], 1.0 - std::tanh(2.0) * std::tanh(2.0));
}
