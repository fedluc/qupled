#ifndef MPI_UTIL_HPP
#define MPI_UTIL_HPP

#include <cassert>
#include <functional>
#include <string>
#include <vector>

/**
 * @brief Utility functions for optional MPI-based parallel calculations.
 *
 * Provides a thin abstraction layer over MPI that compiles to no-op stubs
 * when MPI is disabled (@c USE_MPI not defined).  The @p isUsed flag can be
 * queried at runtime to determine whether MPI support is active.
 */
namespace MPIUtil {

/** @brief True if the library was compiled with MPI support. */
#ifdef USE_MPI
  constexpr bool isUsed = true;
#else
  constexpr bool isUsed = false;
#endif

  /** @brief Initialize the MPI runtime. */
  void init();

  /** @brief Finalize the MPI runtime and release resources. */
  void finalize();

  /**
   * @brief Return true if @c MPI_Init has been called.
   * @return True if MPI is initialized.
   */
  bool isInitialized();

  /**
   * @brief Return the rank of the calling process.
   * @return Rank in [0, @p numberOfRanks() - 1].
   */
  int rank();

  /**
   * @brief Return the total number of MPI processes in @c MPI_COMM_WORLD.
   * @return Process count.
   */
  int numberOfRanks();

  /** @brief Synchronize all MPI processes at a barrier. */
  void barrier();

  /**
   * @brief Return true if the calling process is the root (rank 0).
   * @return True for rank 0.
   */
  bool isRoot();

  /**
   * @brief Return true if only one MPI process is running.
   * @return True if @p numberOfRanks() == 1.
   */
  bool isSingleProcess();

  /**
   * @brief Throw a runtime error with the given message.
   * @param errMsg Human-readable error description.
   */
  void throwError(const std::string &errMsg);

  /** @brief Abort the MPI job immediately. */
  void abort();

  /**
   * @brief Return the wall-clock time in seconds.
   * @return Elapsed wall time.
   */
  double timer();

  /**
   * @brief Check that all ranks hold the same integer value.
   * @param myNumber Local integer value to compare across ranks.
   * @return True if all ranks agree on @p myNumber.
   */
  bool isEqualOnAllRanks(const int &myNumber);

  /**
   * @brief Per-rank loop index range (inclusive start, exclusive end).
   *
   * Used to distribute a loop of @p loopSize iterations across MPI ranks.
   */
  using MPIParallelForData = std::vector<std::pair<int, int>>;

  /**
   * @brief Compute the loop index range for the calling rank.
   * @param loopSize  Total number of loop iterations.
   * @param thisRank  Rank of the calling process.
   * @return Pair {start, end} for this rank's slice.
   */
  std::pair<int, int> getLoopIndexes(const int loopSize, const int thisRank);

  /**
   * @brief Compute loop index ranges for all ranks.
   * @param loopSize Total number of loop iterations.
   * @return Vector of {start, end} pairs, one per rank.
   */
  MPIParallelForData getAllLoopIndexes(const int loopSize);

  /**
   * @brief Execute a parallel for loop across all MPI ranks.
   *
   * Distributes @p loopSize iterations across ranks, invokes @p loopFunc
   * on each rank's slice (optionally parallelized with OpenMP), and returns
   * the per-rank index ranges.
   *
   * @param loopFunc    Function to call with the loop index.
   * @param loopSize    Total number of iterations.
   * @param ompThreads  Number of OpenMP threads per rank.
   * @return Per-rank index ranges (needed for @p gatherLoopData).
   */
  MPIParallelForData parallelFor(const std::function<void(int)> &loopFunc,
                                 const int loopSize,
                                 const int ompThreads);

  /**
   * @brief Gather and synchronize data produced by a parallel for loop.
   *
   * Each rank has computed @p countsPerLoop values for its slice of the loop.
   * After this call @p dataToGather is consistent on all ranks.
   *
   * @param dataToGather  Pointer to the data buffer (size = loopSize ×
   * countsPerLoop).
   * @param loopData      Per-rank index ranges returned by @p parallelFor.
   * @param countsPerLoop Number of values produced per loop iteration.
   */
  void gatherLoopData(double *dataToGather,
                      const MPIParallelForData &loopData,
                      const int countsPerLoop);

} // namespace MPIUtil

#endif
