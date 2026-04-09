#include <gtest/gtest.h>

#include <string>

#include "util/format.hpp"

TEST(FormatTest, FormatsPositionalArgumentsAcrossPlatforms) {
  const std::string out = formatUtil::format("run {}: {:.1f}", 7, 3.25);
  EXPECT_EQ(out, "run 7: 3.2");
}
