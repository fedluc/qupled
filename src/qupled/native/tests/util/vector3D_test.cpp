#include <gtest/gtest.h>

#include <vector>

#include "util/vector3D.hpp"

TEST(Vector3DTest, SizeAccessAndFillOverloadsBehaveAsExpected) {
  Vector3D tensor(2, 2, 3);
  EXPECT_EQ(tensor.size(), 12u);
  EXPECT_EQ(tensor.size(0), 2u);
  EXPECT_EQ(tensor.size(1), 2u);
  EXPECT_EQ(tensor.size(2), 3u);

  tensor.fill(1.0);
  EXPECT_DOUBLE_EQ(tensor(1, 1, 2), 1.0);

  tensor.fill(0, 2.0);
  EXPECT_DOUBLE_EQ(tensor(0, 0, 0), 2.0);
  EXPECT_DOUBLE_EQ(tensor(0, 1, 2), 2.0);

  tensor.fill(1, 0, 3.0);
  EXPECT_DOUBLE_EQ(tensor(1, 0, 0), 3.0);
  EXPECT_DOUBLE_EQ(tensor(1, 0, 2), 3.0);

  tensor.fill(1, 1, std::vector<double>{4.0, 5.0, 6.0});
  EXPECT_DOUBLE_EQ(tensor(1, 1, 0), 4.0);
  EXPECT_DOUBLE_EQ(tensor(1, 1, 2), 6.0);
}

TEST(Vector3DTest, ElementWiseArithmeticAndScalingWork) {
  Vector3D left(1, 2, 2);
  left.fill(1.0);
  Vector3D right(1, 2, 2);
  right.fill(2.0);

  left.sum(right);
  EXPECT_DOUBLE_EQ(left(0, 0, 0), 3.0);

  left.diff(right);
  EXPECT_DOUBLE_EQ(left(0, 1, 1), 1.0);

  left.mult(right);
  EXPECT_DOUBLE_EQ(left(0, 0, 1), 2.0);

  left.div(right);
  EXPECT_DOUBLE_EQ(left(0, 1, 0), 1.0);

  left.mult(0.25);
  EXPECT_DOUBLE_EQ(left(0, 0, 0), 0.25);

  left.linearCombination(right, -2.0);
  EXPECT_DOUBLE_EQ(left(0, 0, 0), -3.75);
}

TEST(Vector3DTest, ResizeResetsDataAndEqualityChecksDimensions) {
  Vector3D a(2, 1, 2);
  a.fill(7.0);
  Vector3D b(2, 1, 2);
  b.fill(7.0);
  EXPECT_TRUE(a == b);

  a.resize(1, 1, 1);
  EXPECT_EQ(a.size(), 1u);
  EXPECT_DOUBLE_EQ(a(0, 0, 0), 0.0);
  EXPECT_FALSE(a == b);
}
