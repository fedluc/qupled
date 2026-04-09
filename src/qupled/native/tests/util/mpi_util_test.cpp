#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "util/mpi_util.hpp"

TEST(MpiUtilCompileFlagTest, IsUsedFlagIsBooleanLike) {
  EXPECT_TRUE(MPIUtil::isUsed || !MPIUtil::isUsed);
}

TEST(MpiUtilRankTest, NumberOfRanksIsAtLeastOne) {
  EXPECT_GE(MPIUtil::numberOfRanks(), 1);
}

TEST(MpiUtilRankTest, RankIsNonNegative) {
  EXPECT_GE(MPIUtil::rank(), 0);
}

TEST(MpiUtilRankTest, RankIsStrictlySmallerThanNumberOfRanks) {
  EXPECT_LT(MPIUtil::rank(), MPIUtil::numberOfRanks());
}

TEST(MpiUtilRoleTest, IsRootMatchesRankEqualsZero) {
  EXPECT_EQ(MPIUtil::isRoot(), MPIUtil::rank() == 0);
}

TEST(MpiUtilRoleTest, IsSingleProcessMatchesWorldSizeEqualsOne) {
  EXPECT_EQ(MPIUtil::isSingleProcess(), MPIUtil::numberOfRanks() == 1);
}

TEST(MpiUtilLifecycleTest, IsInitializedIsQueryable) {
  EXPECT_TRUE(MPIUtil::isInitialized() || !MPIUtil::isInitialized());
}

TEST(MpiUtilLifecycleTest, InitAndFinalizeAreCallableWhenMpiIsDisabled) {
  if constexpr (!MPIUtil::isUsed) {
    EXPECT_NO_THROW(MPIUtil::init());
    EXPECT_NO_THROW(MPIUtil::finalize());
  } else {
    GTEST_SKIP() << "MPI lifecycle is controlled externally in MPI-enabled runs.";
  }
}

TEST(MpiUtilLifecycleTest, AbortIsCallableWhenMpiIsDisabled) {
  if constexpr (!MPIUtil::isUsed) {
    EXPECT_NO_THROW(MPIUtil::abort());
  } else {
    GTEST_SKIP() << "Abort is intentionally not called in MPI-enabled runs.";
  }
}

TEST(MpiUtilSyncTest, BarrierIsCallable) {
  EXPECT_NO_THROW(MPIUtil::barrier());
}

TEST(MpiUtilTimeTest, TimerIsMonotonicAcrossBarrier) {
  const double t0 = MPIUtil::timer();
  MPIUtil::barrier();
  const double t1 = MPIUtil::timer();
  EXPECT_GE(t1, t0);
}

TEST(MpiUtilReductionTest, IsEqualOnAllRanksReturnsTrueForUniformValue) {
  EXPECT_TRUE(MPIUtil::isEqualOnAllRanks(7));
}

TEST(MpiUtilLoopIndexesTest, GetLoopIndexesStartIsWithinBounds) {
  constexpr int loop_size = 10;
  const auto idx = MPIUtil::getLoopIndexes(loop_size, MPIUtil::rank());
  EXPECT_GE(idx.first, 0);
}

TEST(MpiUtilLoopIndexesTest, GetLoopIndexesEndIsWithinBounds) {
  constexpr int loop_size = 10;
  const auto idx = MPIUtil::getLoopIndexes(loop_size, MPIUtil::rank());
  EXPECT_LE(idx.second, loop_size);
}

TEST(MpiUtilLoopIndexesTest, GetLoopIndexesProducesNonNegativeSpan) {
  constexpr int loop_size = 10;
  const auto idx = MPIUtil::getLoopIndexes(loop_size, MPIUtil::rank());
  EXPECT_GE(idx.second - idx.first, 0);
}

TEST(MpiUtilLoopIndexesTest, GetAllLoopIndexesReturnsOneRangePerRank) {
  const auto all_idx = MPIUtil::getAllLoopIndexes(11);
  EXPECT_EQ(all_idx.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));
}

TEST(MpiUtilLoopIndexesTest, GetAllLoopIndexesStartsFromZero) {
  const auto all_idx = MPIUtil::getAllLoopIndexes(11);
  ASSERT_FALSE(all_idx.empty());
  EXPECT_EQ(all_idx.front().first, 0);
}

TEST(MpiUtilLoopIndexesTest, GetAllLoopIndexesEndsAtLoopSize) {
  constexpr int loop_size = 11;
  const auto all_idx = MPIUtil::getAllLoopIndexes(loop_size);
  ASSERT_FALSE(all_idx.empty());
  EXPECT_EQ(all_idx.back().second, loop_size);
}

TEST(MpiUtilLoopIndexesTest, GetAllLoopIndexesProducesContiguousRanges) {
  constexpr int loop_size = 13;
  const auto all_idx = MPIUtil::getAllLoopIndexes(loop_size);
  for (size_t i = 1; i < all_idx.size(); ++i) {
    EXPECT_EQ(all_idx[i - 1].second, all_idx[i].first);
  }
}

TEST(MpiUtilParallelForTest, ParallelForReturnsOneRangePerRank) {
  std::vector<int> hits(8, 0);
  const auto loop_data = MPIUtil::parallelFor(
      [&](const int i) { hits[i] += 1; },
      static_cast<int>(hits.size()),
      1);
  EXPECT_EQ(loop_data.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));
}

TEST(MpiUtilParallelForTest, ParallelForExecutesEachIndexInLocalRangeOnce) {
  std::vector<int> hits(8, 0);
  const auto loop_data = MPIUtil::parallelFor(
      [&](const int i) { hits[i] += 1; },
      static_cast<int>(hits.size()),
      2);

  const auto my_idx = loop_data[MPIUtil::rank()];
  for (int i = my_idx.first; i < my_idx.second; ++i) {
    EXPECT_EQ(hits[i], 1);
  }
}

TEST(MpiUtilParallelForTest, ParallelForDoesNotTouchIndexesOutsideLocalRange) {
  std::vector<int> hits(8, 0);
  const auto loop_data = MPIUtil::parallelFor(
      [&](const int i) { hits[i] += 1; },
      static_cast<int>(hits.size()),
      2);

  const auto my_idx = loop_data[MPIUtil::rank()];
  for (int i = 0; i < my_idx.first; ++i) {
    EXPECT_EQ(hits[i], 0);
  }
  for (int i = my_idx.second; i < static_cast<int>(hits.size()); ++i) {
    EXPECT_EQ(hits[i], 0);
  }
}

TEST(MpiUtilGatherTest, GatherLoopDataAcceptsValidBuffer) {
  std::vector<double> values{1.0, 2.0, 3.0};
  const auto loop_data = MPIUtil::getAllLoopIndexes(static_cast<int>(values.size()));
  EXPECT_NO_THROW(MPIUtil::gatherLoopData(values.data(), loop_data, 1));
}

TEST(MpiUtilGatherTest, GatherLoopDataRejectsNullBufferInSingleProcessMode) {
  if (MPIUtil::isSingleProcess()) {
    const auto loop_data = MPIUtil::getAllLoopIndexes(0);
    EXPECT_THROW(MPIUtil::gatherLoopData(nullptr, loop_data, 1),
                 std::runtime_error);
  } else {
    GTEST_SKIP() << "Null-buffer path aborts under multi-rank MPI runs.";
  }
}

TEST(MpiUtilErrorTest, ThrowErrorThrowsRuntimeErrorInSingleProcessMode) {
  if (MPIUtil::isSingleProcess()) {
    EXPECT_THROW(MPIUtil::throwError("test failure"), std::runtime_error);
  } else {
    GTEST_SKIP() << "throwError aborts under multi-rank MPI runs.";
  }
}
