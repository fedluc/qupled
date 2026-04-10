#include <gtest/gtest.h>

#include "vs/grid_point.hpp"

TEST(GridPointTest, RsDownThetaDownMapsToIndexZero) {
  EXPECT_EQ(GridPoints::RS_DOWN_THETA_DOWN.toIndex(), 0u);
}

TEST(GridPointTest, RsThetaDownMapsToIndexOne) {
  EXPECT_EQ(GridPoints::RS_THETA_DOWN.toIndex(), 1u);
}

TEST(GridPointTest, RsUpThetaDownMapsToIndexTwo) {
  EXPECT_EQ(GridPoints::RS_UP_THETA_DOWN.toIndex(), 2u);
}

TEST(GridPointTest, RsDownThetaMapsToIndexThree) {
  EXPECT_EQ(GridPoints::RS_DOWN_THETA.toIndex(), 3u);
}

TEST(GridPointTest, CenterMapsToIndexFour) {
  EXPECT_EQ(GridPoints::CENTER.toIndex(), 4u);
}

TEST(GridPointTest, RsUpThetaMapsToIndexFive) {
  EXPECT_EQ(GridPoints::RS_UP_THETA.toIndex(), 5u);
}

TEST(GridPointTest, RsDownThetaUpMapsToIndexSix) {
  EXPECT_EQ(GridPoints::RS_DOWN_THETA_UP.toIndex(), 6u);
}

TEST(GridPointTest, RsThetaUpMapsToIndexSeven) {
  EXPECT_EQ(GridPoints::RS_THETA_UP.toIndex(), 7u);
}

TEST(GridPointTest, RsUpThetaUpMapsToIndexEight) {
  EXPECT_EQ(GridPoints::RS_UP_THETA_UP.toIndex(), 8u);
}

TEST(GridPointTest, UsesThetaAsOuterIndexAndRsAsInnerIndex) {
  const GridPoint custom{GridPoint::Rs::UP, GridPoint::Theta::DOWN};
  EXPECT_EQ(custom.toIndex(), 2u);
}
