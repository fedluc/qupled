#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "util/mpi_util.hpp"

TEST(MpiUtilTest, ExposesACompileTimeMpiUsageFlag) {
  EXPECT_TRUE(MPIUtil::isUsed || !MPIUtil::isUsed);
}

TEST(MpiUtilTest, InitAndFinalizeAreNoOpsWhenMpiIsDisabled) {
  if constexpr (!MPIUtil::isUsed) {
    EXPECT_NO_THROW(MPIUtil::init());
    EXPECT_NO_THROW(MPIUtil::finalize());
  } else {
    GTEST_SKIP() << "Lifecycle is managed externally when MPI is enabled.";
  }
}

TEST(MpiUtilTest, IsInitializedReportsTrueWhenMpiIsDisabled) {
  if constexpr (!MPIUtil::isUsed) {
    EXPECT_TRUE(MPIUtil::isInitialized());
  } else {
    SUCCEED() << "Initialization status depends on MPI runtime lifecycle.";
  }
}

TEST(MpiUtilTest, RankAndWorldSizeStayInValidRange) {
  const int rank = MPIUtil::rank();
  const int n_ranks = MPIUtil::numberOfRanks();

  EXPECT_GE(n_ranks, 1);
  EXPECT_GE(rank, 0);
  EXPECT_LT(rank, n_ranks);
}

TEST(MpiUtilTest, RootFlagMatchesRankZero) {
  EXPECT_EQ(MPIUtil::isRoot(), MPIUtil::rank() == 0);
}

TEST(MpiUtilTest, SingleProcessFlagMatchesWorldSize) {
  EXPECT_EQ(MPIUtil::isSingleProcess(), MPIUtil::numberOfRanks() == 1);
}

TEST(MpiUtilTest, BarrierCanBeInvokedSafely) {
  EXPECT_NO_THROW(MPIUtil::barrier());
}

TEST(MpiUtilTest, TimerDoesNotGoBackwards) {
  const double t0 = MPIUtil::timer();
  MPIUtil::barrier();
  const double t1 = MPIUtil::timer();
  EXPECT_GE(t1, t0);
}

TEST(MpiUtilTest, IsEqualOnAllRanksAcceptsUniformValue) {
  EXPECT_TRUE(MPIUtil::isEqualOnAllRanks(42));
}

TEST(MpiUtilTest, GetLoopIndexesReturnsLocalBoundsInsideLoopRange) {
  constexpr int loop_size = 10;
  const auto idx = MPIUtil::getLoopIndexes(loop_size, MPIUtil::rank());

  EXPECT_GE(idx.first, 0);
  EXPECT_LE(idx.first, loop_size);
  EXPECT_GE(idx.second, idx.first);
  EXPECT_LE(idx.second, loop_size);
}

TEST(MpiUtilTest, GetAllLoopIndexesCoversWholeLoopWithoutGaps) {
  constexpr int loop_size = 13;
  const auto all_idx = MPIUtil::getAllLoopIndexes(loop_size);

  ASSERT_EQ(all_idx.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));

  int next_start = 0;
  for (const auto &idx : all_idx) {
    EXPECT_EQ(idx.first, next_start);
    EXPECT_GE(idx.second, idx.first);
    EXPECT_LE(idx.second, loop_size);
    next_start = idx.second;
  }
  EXPECT_EQ(next_start, loop_size);
}

TEST(MpiUtilTest, ParallelForExecutesOnlyThisRanksAssignedSlice) {
  std::vector<int> hits(8, 0);
  const auto loop_data = MPIUtil::parallelFor(
      [&](const int i) { hits[i] += 1; },
      static_cast<int>(hits.size()),
      2);

  ASSERT_EQ(loop_data.size(), static_cast<size_t>(MPIUtil::numberOfRanks()));
  const auto my_idx = loop_data[MPIUtil::rank()];
  for (int i = 0; i < static_cast<int>(hits.size()); ++i) {
    if (i >= my_idx.first && i < my_idx.second) {
      EXPECT_EQ(hits[i], 1);
    } else {
      EXPECT_EQ(hits[i], 0);
    }
  }
}

TEST(MpiUtilTest, GatherLoopDataAcceptsValidBuffer) {
  std::vector<double> values{1.0, 2.0, 3.0};
  const auto loop_data = MPIUtil::getAllLoopIndexes(static_cast<int>(values.size()));
  EXPECT_NO_THROW(MPIUtil::gatherLoopData(values.data(), loop_data, 1));
}

TEST(MpiUtilTest, GatherLoopDataRejectsNullBufferInSingleProcessMode) {
  if (MPIUtil::isSingleProcess()) {
    const auto loop_data = MPIUtil::getAllLoopIndexes(0);
    EXPECT_THROW(MPIUtil::gatherLoopData(nullptr, loop_data, 1),
                 std::runtime_error);
  } else {
    GTEST_SKIP() << "Failure path aborts MPI jobs when multiple ranks are used.";
  }
}

TEST(MpiUtilTest, ThrowErrorRaisesRuntimeErrorInSingleProcessMode) {
  if (MPIUtil::isSingleProcess()) {
    EXPECT_THROW(MPIUtil::throwError("test failure"), std::runtime_error);
  } else {
    GTEST_SKIP() << "Failure path aborts MPI jobs when multiple ranks are used.";
  }
}

TEST(MpiUtilTest, AbortIsNoOpWhenMpiIsDisabled) {
  if constexpr (!MPIUtil::isUsed) {
    EXPECT_NO_THROW(MPIUtil::abort());
  } else {
    GTEST_SKIP() << "Abort is intentionally not invoked in MPI-enabled tests.";
  }
}
