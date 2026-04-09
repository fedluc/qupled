#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "util/vector_util.hpp"

TEST(VectorUtilTest, SumReturnsElementwiseResult) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};

  const auto out = vecUtil::sum(v1, v2);

  ASSERT_EQ(out.size(), 3u);
  EXPECT_DOUBLE_EQ(out[0], 5.0);
  EXPECT_DOUBLE_EQ(out[1], 7.0);
  EXPECT_DOUBLE_EQ(out[2], 9.0);
}

TEST(VectorUtilTest, DiffReturnsElementwiseResult) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};

  const auto out = vecUtil::diff(v1, v2);

  ASSERT_EQ(out.size(), 3u);
  EXPECT_DOUBLE_EQ(out[0], 3.0);
  EXPECT_DOUBLE_EQ(out[1], 3.0);
  EXPECT_DOUBLE_EQ(out[2], 3.0);
}

TEST(VectorUtilTest, MultReturnsElementwiseProduct) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};

  const auto out = vecUtil::mult(v1, v2);

  ASSERT_EQ(out.size(), 3u);
  EXPECT_DOUBLE_EQ(out[0], 4.0);
  EXPECT_DOUBLE_EQ(out[1], 10.0);
  EXPECT_DOUBLE_EQ(out[2], 18.0);
}

TEST(VectorUtilTest, DivReturnsElementwiseQuotient) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};

  const auto out = vecUtil::div(v1, v2);

  ASSERT_EQ(out.size(), 3u);
  EXPECT_DOUBLE_EQ(out[0], 4.0);
  EXPECT_DOUBLE_EQ(out[1], 2.5);
  EXPECT_DOUBLE_EQ(out[2], 2.0);
}

TEST(VectorUtilTest, MultWithScalarScalesAllEntries) {
  const std::vector<double> v1{1.0, 2.0, 3.0};

  const auto scaled = vecUtil::mult(v1, 2.5);

  ASSERT_EQ(scaled.size(), 3u);
  EXPECT_DOUBLE_EQ(scaled[0], 2.5);
  EXPECT_DOUBLE_EQ(scaled[1], 5.0);
  EXPECT_DOUBLE_EQ(scaled[2], 7.5);
}

TEST(VectorUtilTest, LinearCombinationReturnsWeightedSum) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{2.0, 4.0, 8.0};

  const auto linear = vecUtil::linearCombination(v1, 2.0, v2, -0.5);

  ASSERT_EQ(linear.size(), 3u);
  EXPECT_DOUBLE_EQ(linear[0], 1.0);
  EXPECT_DOUBLE_EQ(linear[1], 2.0);
  EXPECT_DOUBLE_EQ(linear[2], 2.0);
}

TEST(VectorUtilTest, RmsReturnsRootSumSquaredErrorWhenNotNormalized) {
  const std::vector<double> values{1.0, -2.0, 3.5};
  const std::vector<double> zeros{0.0, 0.0, 0.0};

  const double out = vecUtil::rms(values, zeros, false);
  EXPECT_DOUBLE_EQ(out, std::sqrt(1.0 + 4.0 + 12.25));
}

TEST(VectorUtilTest, RmsReturnsRootMeanSquaredErrorWhenNormalized) {
  const std::vector<double> values{1.0, -2.0, 3.5};
  const std::vector<double> zeros{0.0, 0.0, 0.0};

  const double out = vecUtil::rms(values, zeros, true);
  EXPECT_DOUBLE_EQ(out, std::sqrt((1.0 + 4.0 + 12.25) / 3.0));
}

TEST(VectorUtilTest, FillOverwritesEachElement) {
  std::vector<double> values{1.0, -2.0, 3.5};

  vecUtil::fill(values, 7.0);

  ASSERT_EQ(values.size(), 3u);
  EXPECT_DOUBLE_EQ(values[0], 7.0);
  EXPECT_DOUBLE_EQ(values[1], 7.0);
  EXPECT_DOUBLE_EQ(values[2], 7.0);
}

#ifndef NDEBUG
TEST(VectorUtilTest, SumDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::sum(v1, v2), "");
}

TEST(VectorUtilTest, LinearCombinationDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::linearCombination(v1, 1.0, v2, 1.0), "");
}
#endif
