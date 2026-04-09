#include <gtest/gtest.h>

#include "vs/grid_point.hpp"

TEST(GridPointTest, NamedPointsMapToStableFlatIndexes) {
  EXPECT_EQ(GridPoints::RS_DOWN_THETA_DOWN.toIndex(), 0u);
  EXPECT_EQ(GridPoints::RS_THETA_DOWN.toIndex(), 1u);
  EXPECT_EQ(GridPoints::RS_UP_THETA_DOWN.toIndex(), 2u);
  EXPECT_EQ(GridPoints::RS_DOWN_THETA.toIndex(), 3u);
  EXPECT_EQ(GridPoints::CENTER.toIndex(), 4u);
  EXPECT_EQ(GridPoints::RS_UP_THETA.toIndex(), 5u);
  EXPECT_EQ(GridPoints::RS_DOWN_THETA_UP.toIndex(), 6u);
  EXPECT_EQ(GridPoints::RS_THETA_UP.toIndex(), 7u);
  EXPECT_EQ(GridPoints::RS_UP_THETA_UP.toIndex(), 8u);
}
