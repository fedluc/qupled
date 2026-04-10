#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "util/vector_util.hpp"

TEST(VectorSumTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::sum(v1, v2);
  EXPECT_DOUBLE_EQ(out[0], 5.0);
}

TEST(VectorSumTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::sum(v1, v2);
  EXPECT_DOUBLE_EQ(out[1], 7.0);
}

TEST(VectorSumTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::sum(v1, v2);
  EXPECT_DOUBLE_EQ(out[2], 9.0);
}

TEST(VectorSumTest, ReturnsOutputWithSameSizeAsInputs) {
  const std::vector<double> v1{1.0, 2.0, 3.0, 4.0};
  const std::vector<double> v2{4.0, 3.0, 2.0, 1.0};
  const auto out = vecUtil::sum(v1, v2);
  EXPECT_EQ(out.size(), v1.size());
}

TEST(VectorDiffTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::diff(v1, v2);
  EXPECT_DOUBLE_EQ(out[0], 3.0);
}

TEST(VectorDiffTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::diff(v1, v2);
  EXPECT_DOUBLE_EQ(out[1], 3.0);
}

TEST(VectorDiffTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::diff(v1, v2);
  EXPECT_DOUBLE_EQ(out[2], 3.0);
}

TEST(VectorMultTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::mult(v1, v2);
  EXPECT_DOUBLE_EQ(out[0], 4.0);
}

TEST(VectorMultTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::mult(v1, v2);
  EXPECT_DOUBLE_EQ(out[1], 10.0);
}

TEST(VectorMultTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{4.0, 5.0, 6.0};
  const auto out = vecUtil::mult(v1, v2);
  EXPECT_DOUBLE_EQ(out[2], 18.0);
}

TEST(VectorDivTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::div(v1, v2);
  EXPECT_DOUBLE_EQ(out[0], 4.0);
}

TEST(VectorDivTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::div(v1, v2);
  EXPECT_DOUBLE_EQ(out[1], 2.5);
}

TEST(VectorDivTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{4.0, 5.0, 6.0};
  const std::vector<double> v2{1.0, 2.0, 3.0};
  const auto out = vecUtil::div(v1, v2);
  EXPECT_DOUBLE_EQ(out[2], 2.0);
}

TEST(VectorScalarMultTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const auto scaled = vecUtil::mult(v1, 2.5);
  EXPECT_DOUBLE_EQ(scaled[0], 2.5);
}

TEST(VectorScalarMultTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const auto scaled = vecUtil::mult(v1, 2.5);
  EXPECT_DOUBLE_EQ(scaled[1], 5.0);
}

TEST(VectorScalarMultTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const auto scaled = vecUtil::mult(v1, 2.5);
  EXPECT_DOUBLE_EQ(scaled[2], 7.5);
}

TEST(VectorLinearCombinationTest, ReturnsExpectedFirstEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{2.0, 4.0, 8.0};
  const auto linear = vecUtil::linearCombination(v1, 2.0, v2, -0.5);
  EXPECT_DOUBLE_EQ(linear[0], 1.0);
}

TEST(VectorLinearCombinationTest, ReturnsExpectedSecondEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{2.0, 4.0, 8.0};
  const auto linear = vecUtil::linearCombination(v1, 2.0, v2, -0.5);
  EXPECT_DOUBLE_EQ(linear[1], 2.0);
}

TEST(VectorLinearCombinationTest, ReturnsExpectedThirdEntry) {
  const std::vector<double> v1{1.0, 2.0, 3.0};
  const std::vector<double> v2{2.0, 4.0, 8.0};
  const auto linear = vecUtil::linearCombination(v1, 2.0, v2, -0.5);
  EXPECT_DOUBLE_EQ(linear[2], 2.0);
}

TEST(VectorRmsTest, ReturnsRootSumSquaredErrorWhenNormalizeIsFalse) {
  const std::vector<double> values{1.0, -2.0, 3.5};
  const std::vector<double> zeros{0.0, 0.0, 0.0};
  const double out = vecUtil::rms(values, zeros, false);
  EXPECT_DOUBLE_EQ(out, std::sqrt(1.0 + 4.0 + 12.25));
}

TEST(VectorRmsTest, ReturnsRootMeanSquaredErrorWhenNormalizeIsTrue) {
  const std::vector<double> values{1.0, -2.0, 3.5};
  const std::vector<double> zeros{0.0, 0.0, 0.0};
  const double out = vecUtil::rms(values, zeros, true);
  EXPECT_DOUBLE_EQ(out, std::sqrt((1.0 + 4.0 + 12.25) / 3.0));
}

TEST(VectorFillTest, OverwritesFirstEntry) {
  std::vector<double> values{1.0, -2.0, 3.5};
  vecUtil::fill(values, 7.0);
  EXPECT_DOUBLE_EQ(values[0], 7.0);
}

TEST(VectorFillTest, OverwritesSecondEntry) {
  std::vector<double> values{1.0, -2.0, 3.5};
  vecUtil::fill(values, 7.0);
  EXPECT_DOUBLE_EQ(values[1], 7.0);
}

TEST(VectorFillTest, OverwritesThirdEntry) {
  std::vector<double> values{1.0, -2.0, 3.5};
  vecUtil::fill(values, 7.0);
  EXPECT_DOUBLE_EQ(values[2], 7.0);
}

#ifndef NDEBUG
TEST(VectorSizeMismatchTest, SumDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::sum(v1, v2), "");
}

TEST(VectorSizeMismatchTest, DiffDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::diff(v1, v2), "");
}

TEST(VectorSizeMismatchTest, MultDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::mult(v1, v2), "");
}

TEST(VectorSizeMismatchTest, DivDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::div(v1, v2), "");
}

TEST(VectorSizeMismatchTest, LinearCombinationDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::linearCombination(v1, 1.0, v2, 1.0), "");
}

TEST(VectorSizeMismatchTest, RmsDiesWhenVectorSizesDiffer) {
  const std::vector<double> v1{1.0, 2.0};
  const std::vector<double> v2{1.0};
  EXPECT_DEATH((void)vecUtil::rms(v1, v2, false), "");
}
#endif
