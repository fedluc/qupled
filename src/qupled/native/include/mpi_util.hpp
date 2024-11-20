#ifndef MPI_UTIL_HPP
#define MPI_UTIL_HPP

#include <cassert>
#include <functional>
#include <string>
#include <vector>

// -------------------------------------------------------------------
// Utility functions to handle parallel calculations with MPI
// -------------------------------------------------------------------

namespace MPIUtil {

  // Initialize MPI
  void init();

  // Finalize MPI
  void finalize();

  // Check if MPI initialized
  bool isInitialized();

  // Get rank of MPI process
  int rank();

  // Get total number of MPI processes
  int numberOfRanks();

  // Set an MPI Barrier
  void barrier();

  // Check if the process is the root process
  bool isRoot();

  // Check if only one rank is used
  bool isSingleProcess();

  // Throw error with description given in errMsg
  void throwError(const std::string &errMsg);

  // Abort MPI
  void abort();

  // Get wall time
  double timer();

  // Check that a number is the same on all ranks
  bool isEqualOnAllRanks(const int &myNumber);

  // Data structure to track how loop indexes are distributed
  using MPIParallelForData = std::vector<std::pair<int, int>>;

  // Get start and finish index for parallel for loop on one rank
  std::pair<int, int> getLoopIndexes(const int loopSize, const int thisRank);

  // Get start and finish index for parallel for loop on all ranks
  MPIParallelForData getAllLoopIndexes(const int loopSize);

  // Wrapper for parallel for loop
  MPIParallelForData parallelFor(const std::function<void(int)> &loopFunc,
                                 const int loopSize,
                                 const int ompThreads);

  // Synchronize data from a parallel for loop among all ranks
  void gatherLoopData(double *dataToGather,
                      const MPIParallelForData &loopData,
                      const int countsPerLoop);

} // namespace MPIUtil

#endif
