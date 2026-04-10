#include <gtest/gtest.h>

#include <stdexcept>

#include "util/dimensions_util.hpp"

namespace {
  class FakeDimensionsHandler : public dimensionsUtil::DimensionsHandler {
  public:

    int compute2d_calls = 0;
    int compute3d_calls = 0;

  private:

    void compute2D() override { ++compute2d_calls; }
    void compute3D() override { ++compute3d_calls; }
  };
} // namespace

TEST(DimensionsUtilTest, DispatchesD2To2DImplementation) {
  FakeDimensionsHandler handler;
  handler.compute(dimensionsUtil::Dimension::D2);
  EXPECT_EQ(handler.compute2d_calls, 1);
}

TEST(DimensionsUtilTest, D2DispatchDoesNotCall3DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::D2);
  EXPECT_EQ(handler.compute3d_calls, 0);
}

TEST(DimensionsUtilTest, DispatchesD3To3DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::D3);
  EXPECT_EQ(handler.compute3d_calls, 1);
}

TEST(DimensionsUtilTest, D3DispatchDoesNotCall2DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::D3);
  EXPECT_EQ(handler.compute2d_calls, 0);
}

TEST(DimensionsUtilTest, DispatchesDefaultDimensionTo3DImplementation) {
  FakeDimensionsHandler handler;
  handler.compute(dimensionsUtil::Dimension::Default);
  EXPECT_EQ(handler.compute3d_calls, 1);
}

TEST(DimensionsUtilTest, DefaultDimensionDispatchDoesNotCall2DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::Default);
  EXPECT_EQ(handler.compute2d_calls, 0);
}

TEST(DimensionsUtilTest, ThrowsOnUnsupportedDimension) {
  FakeDimensionsHandler handler;
  const auto invalid = static_cast<dimensionsUtil::Dimension>(12345);
  EXPECT_THROW(handler.compute(invalid), std::runtime_error);
}
