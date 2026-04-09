#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "util/mpi_util.hpp"

TEST(MpiUtilTest, ReportsConsistentSingleProcessBasics) {
  EXPECT_TRUE(MPIUtil::isInitialized());
  EXPECT_EQ(MPIUtil::rank(), 0);
  EXPECT_EQ(MPIUtil::numberOfRanks(), 1);
  EXPECT_TRUE(MPIUtil::isRoot());
  EXPECT_TRUE(MPIUtil::isSingleProcess());
}

TEST(MpiUtilTest, ThrowErrorRaisesRuntimeErrorInSingleProcessMode) {
  EXPECT_THROW(MPIUtil::throwError("test failure"), std::runtime_error);
}

TEST(MpiUtilTest, LoopIndexHelpersCoverEntireRange) {
  const auto idx = MPIUtil::getLoopIndexes(10, 0);
  EXPECT_EQ(idx.first, 0);
  EXPECT_EQ(idx.second, 10);

  const auto all_idx = MPIUtil::getAllLoopIndexes(10);
  ASSERT_EQ(all_idx.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));
  EXPECT_EQ(all_idx[0].first, 0);
  EXPECT_LE(all_idx[0].second, 10);
}

TEST(MpiUtilTest, ParallelForExecutesAllIndexesAndReturnsRanges) {
  std::vector<int> hits(8, 0);
  const auto loop_data = MPIUtil::parallelFor(
      [&](const int i) { hits[i] += 1; },
      static_cast<int>(hits.size()),
      1);

  ASSERT_EQ(loop_data.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));
  EXPECT_EQ(loop_data[0].first, 0);
  EXPECT_LE(loop_data[0].second, 8);
  for (const int hit : hits) {
    EXPECT_EQ(hit, 1);
  }
}

TEST(MpiUtilTest, GatherLoopDataAcceptsValidBuffer) {
  std::vector<double> values{1.0, 2.0, 3.0};
  const auto loop_data = MPIUtil::getAllLoopIndexes(static_cast<int>(values.size()));
  EXPECT_NO_THROW(MPIUtil::gatherLoopData(values.data(), loop_data, 1));
}
