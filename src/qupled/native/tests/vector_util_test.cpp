#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "util/vector_util.hpp"

TEST(VectorUtilTest, ElementWiseOperationsProduceExpectedValues) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};

  const auto summed = vecUtil::sum(v1, v2);
  const auto diffed = vecUtil::diff(v2, v1);
  const auto multiplied = vecUtil::mult(v1, v2);
  const auto divided = vecUtil::div(v2, v1);

  ASSERT_EQ(summed.size(), 3u);
  ASSERT_EQ(diffed.size(), 3u);
  ASSERT_EQ(multiplied.size(), 3u);
  ASSERT_EQ(divided.size(), 3u);

  EXPECT_DOUBLE_EQ(summed[0], 5.0);
  EXPECT_DOUBLE_EQ(summed[1], 7.0);
  EXPECT_DOUBLE_EQ(summed[2], 9.0);
  EXPECT_DOUBLE_EQ(diffed[0], 3.0);
  EXPECT_DOUBLE_EQ(diffed[1], 3.0);
  EXPECT_DOUBLE_EQ(diffed[2], 3.0);
  EXPECT_DOUBLE_EQ(multiplied[0], 4.0);
  EXPECT_DOUBLE_EQ(multiplied[1], 10.0);
  EXPECT_DOUBLE_EQ(multiplied[2], 18.0);
  EXPECT_DOUBLE_EQ(divided[0], 4.0);
  EXPECT_DOUBLE_EQ(divided[1], 2.5);
  EXPECT_DOUBLE_EQ(divided[2], 2.0);
}

TEST(VectorUtilTest, LinearCombinationAndScalarMultiplyWork) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{2.0, 4.0, 8.0};

  const auto scaled = vecUtil::mult(v1, 2.5);
  const auto linear = vecUtil::linearCombination(v1, 2.0, v2, -0.5);

  ASSERT_EQ(scaled.size(), 3u);
  ASSERT_EQ(linear.size(), 3u);

  EXPECT_DOUBLE_EQ(scaled[0], 2.5);
  EXPECT_DOUBLE_EQ(scaled[1], 5.0);
  EXPECT_DOUBLE_EQ(scaled[2], 7.5);
  EXPECT_DOUBLE_EQ(linear[0], 1.0);
  EXPECT_DOUBLE_EQ(linear[1], 2.0);
  EXPECT_DOUBLE_EQ(linear[2], 2.0);
}

TEST(VectorUtilTest, RmsAndFillBehaveAsExpected) {
  std::vector<double> values{1.0, -2.0, 3.5};
  const std::vector<double> zeros{0.0, 0.0, 0.0};

  const double rms = vecUtil::rms(values, zeros, false);
  EXPECT_DOUBLE_EQ(rms, std::sqrt(1.0 + 4.0 + 12.25));

  vecUtil::fill(values, 7.0);
  ASSERT_EQ(values.size(), 3u);
  EXPECT_DOUBLE_EQ(values[0], 7.0);
  EXPECT_DOUBLE_EQ(values[1], 7.0);
  EXPECT_DOUBLE_EQ(values[2], 7.0);
}
