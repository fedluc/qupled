#include <gtest/gtest.h>

#include <vector>

#include "util/vector2D.hpp"

TEST(Vector2DTest, ConstructsFromNestedAndFlatVectors) {
  const Vector2D nested(
      std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}});
  EXPECT_EQ(nested.size(), 4u);
  EXPECT_EQ(nested.size(0), 2u);
  EXPECT_EQ(nested.size(1), 2u);
  EXPECT_DOUBLE_EQ(nested(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(nested(1, 0), 3.0);

  const Vector2D flat(std::vector<double>{5.0, 6.0, 7.0});
  EXPECT_EQ(flat.size(0), 3u);
  EXPECT_EQ(flat.size(1), 1u);
  EXPECT_DOUBLE_EQ(flat(2, 0), 7.0);
}

TEST(Vector2DTest, RowSpanAndFillOperationsMutateExpectedEntries) {
  Vector2D matrix(2, 3);
  matrix.fill(1.0);

  auto row0 = matrix[0];
  row0[1] = 9.0;
  EXPECT_DOUBLE_EQ(matrix(0, 1), 9.0);

  matrix.fill(1, -2.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), -2.0);
  EXPECT_DOUBLE_EQ(matrix(1, 2), -2.0);

  matrix.fill(0, std::vector<double>{7.0, 8.0, 9.0});
  EXPECT_DOUBLE_EQ(matrix(0, 0), 7.0);
  EXPECT_DOUBLE_EQ(matrix(0, 2), 9.0);
}

TEST(Vector2DTest, ElementWiseArithmeticAndLinearCombinationWork) {
  Vector2D left(std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}});
  const Vector2D right(
      std::vector<std::vector<double>>{{5.0, 6.0}, {7.0, 8.0}});

  left.sum(right);
  EXPECT_DOUBLE_EQ(left(0, 0), 6.0);
  EXPECT_DOUBLE_EQ(left(1, 1), 12.0);

  left.diff(right);
  EXPECT_DOUBLE_EQ(left(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(left(1, 1), 4.0);

  left.mult(right);
  EXPECT_DOUBLE_EQ(left(0, 0), 5.0);
  EXPECT_DOUBLE_EQ(left(1, 1), 32.0);

  left.div(right);
  EXPECT_DOUBLE_EQ(left(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(left(1, 1), 4.0);

  left.mult(0.5);
  EXPECT_DOUBLE_EQ(left(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(left(1, 1), 2.0);

  left.linearCombination(right, 2.0);
  EXPECT_DOUBLE_EQ(left(0, 0), 10.5);
  EXPECT_DOUBLE_EQ(left(1, 1), 18.0);
}

TEST(Vector2DTest, ResizeReinitializesStorageAndEqualityChecksShape) {
  Vector2D a(2, 2);
  a.fill(3.0);
  Vector2D b(2, 2);
  b.fill(3.0);
  EXPECT_TRUE(a == b);

  a.resize(1, 3);
  EXPECT_EQ(a.size(), 3u);
  EXPECT_EQ(a.size(0), 1u);
  EXPECT_EQ(a.size(1), 3u);
  EXPECT_DOUBLE_EQ(a(0, 0), 0.0);
  EXPECT_FALSE(a == b);
}
