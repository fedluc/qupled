#include "mpi_util.hpp"
#ifdef USE_MPI
#include <mpi.h>
#define OMPI_SKIP_MPICXX 1 // Disable MPI-C++ bindings
#endif
#include <iostream>
#include <numeric>
#include <omp.h>

using namespace std;

namespace MPIUtil {

#ifdef USE_MPI
  const MPI_Comm MPICommunicator = MPI_COMM_WORLD;
#endif

  void init() {
#ifdef USE_MPI
    MPI_Init(nullptr, nullptr);
#endif
  }

  void finalize() {
#ifdef USE_MPI
    MPI_Finalize();
#endif
  }

  bool isInitialized() {
#ifdef USE_MPI
    int isMPIInit;
    MPI_Initialized(&isMPIInit);
    return isMPIInit == 1;
#endif
    return true;
  }

  int rank() {
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPICommunicator, &rank);
    return rank;
#endif
    return 0;
  }

  int numberOfRanks() {
#ifdef USE_MPI
    int numRanks;
    MPI_Comm_size(MPICommunicator, &numRanks);
    return numRanks;
#endif
    return 1;
  }

  void barrier() {
#ifdef USE_MPI
    MPI_Barrier(MPICommunicator);
#endif
  }

  bool isRoot() {
#ifdef USE_MPI
    return rank() == 0;
#endif
    return true;
  }

  bool isSingleProcess() {
#ifdef USE_MPI
    return numberOfRanks() == 1;
#endif
    return true;
  }

  void throwError(const string &errMsg) {
    if (isSingleProcess()) {
      // Throw a catchable error if only one process is used
      throw runtime_error(errMsg);
    }
#ifdef USE_MPI
    // Abort MPI if more than one process is running
    cerr << errMsg << endl;
    abort();
#endif
  }

  void abort() {
#ifdef USE_MPI
    MPI_Abort(MPICommunicator, 1);
#endif
  }

  double timer() {
#ifdef USE_MPI
    return MPI_Wtime();
#endif
    return omp_get_wtime();
  }

  bool isEqualOnAllRanks(const int &myNumber) {
#ifdef USE_MPI
    int globalMininumNumber;
    MPI_Allreduce(
        &myNumber, &globalMininumNumber, 1, MPI_INT, MPI_MIN, MPICommunicator);
    return myNumber == globalMininumNumber;
#endif
    (void)myNumber;
    return true;
  }

  pair<int, int> getLoopIndexes(const int loopSize, const int thisRank) {
    pair<int, int> idx = {0, loopSize};
    const int nRanks = numberOfRanks();
    if (nRanks == 1) { return idx; }
    int localSize = loopSize / nRanks;
    int remainder = loopSize % nRanks;
    idx.first = thisRank * localSize + std::min(thisRank, remainder);
    idx.second = idx.first + localSize + (thisRank < remainder ? 1 : 0);
    idx.second = std::min(idx.second, loopSize);
    return idx;
  }

  MPIParallelForData getAllLoopIndexes(const int loopSize) {
    std::vector<pair<int, int>> out;
    for (int i = 0; i < numberOfRanks(); ++i) {
      out.push_back(getLoopIndexes(loopSize, i));
    }
    return out;
  }

  MPIParallelForData parallelFor(const function<void(int)> &loopFunc,
                                 const int loopSize,
                                 const int ompThreads) {
    MPIParallelForData allIdx = getAllLoopIndexes(loopSize);
    const auto &thisIdx = allIdx[rank()];
    const bool useOMP = ompThreads > 1;
#pragma omp parallel for num_threads(ompThreads) if (useOMP)
    for (int i = thisIdx.first; i < thisIdx.second; ++i) {
      loopFunc(i);
    }
    return allIdx;
  }

  void gatherLoopData(double *dataToGather,
                      const MPIParallelForData &loopData,
                      const int countsPerLoop) {
#ifdef USE_MPI
    std::vector<int> recieverCounts;
    for (const auto &i : loopData) {
      const int loopSpan = i.second - i.first;
      recieverCounts.push_back(loopSpan * countsPerLoop);
    }
    std::vector<int> displacements(recieverCounts.size(), 0);
    std::partial_sum(recieverCounts.begin(),
                     recieverCounts.end() - 1,
                     displacements.begin() + 1,
                     plus<double>());
    MPI_Allgatherv(MPI_IN_PLACE,
                   0,
                   MPI_DATATYPE_NULL,
                   dataToGather,
                   recieverCounts.data(),
                   displacements.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);
#endif
    if (!dataToGather) { throwError(""); }
    (void)loopData;
    (void)countsPerLoop;
  }

} // namespace MPIUtil
