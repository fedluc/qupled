#include <gtest/gtest.h>

#include <stdexcept>

#include "util/dimensions_util.hpp"

namespace {
class FakeDimensionsHandler : public dimensionsUtil::DimensionsHandler {
public:
  int d2_calls = 0;
  int d3_calls = 0;

private:
  void compute2D() override { ++d2_calls; }
  void compute3D() override { ++d3_calls; }
};
} // namespace

TEST(DimensionsUtilTest, DispatchesD2To2DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::D2);

  EXPECT_EQ(handler.d2_calls, 1);
  EXPECT_EQ(handler.d3_calls, 0);
}

TEST(DimensionsUtilTest, DispatchesD3To3DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::D3);

  EXPECT_EQ(handler.d2_calls, 0);
  EXPECT_EQ(handler.d3_calls, 1);
}

TEST(DimensionsUtilTest, DispatchesDefaultDimensionTo3DImplementation) {
  FakeDimensionsHandler handler;

  handler.compute(dimensionsUtil::Dimension::Default);

  EXPECT_EQ(handler.d2_calls, 0);
  EXPECT_EQ(handler.d3_calls, 1);
}

TEST(DimensionsUtilTest, ThrowsOnUnsupportedDimension) {
  FakeDimensionsHandler handler;
  const auto invalid = static_cast<dimensionsUtil::Dimension>(12345);
  EXPECT_THROW(handler.compute(invalid), std::runtime_error);
}
