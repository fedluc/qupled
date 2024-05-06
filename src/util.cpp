#ifdef USE_MPI
#include <mpi.h>
#define OMPI_SKIP_MPICXX 1 // Disable MPI-C++ bindings
#endif

#include <cassert>
#include <numeric>
#include <omp.h>
#include "numerics.hpp"
#include "util.hpp"

using namespace std;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;

namespace thermoUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods
  // -----------------------------------------------------------------

  double computeInternalEnergy(const vector<double> &wvg,
                               const vector<double> &ssf,
                               const double &coupling) {
    const Interpolator1D itp(wvg, ssf);
    Integrator1D itg(1.0e-6);
    const InternalEnergy uInt(coupling, wvg.front(), wvg.back(), itp, itg);
    return uInt.get();
  }

  double computeFreeEnergy(const vector<double> &grid,
                           const vector<double> &rsu,
                           const double &coupling) {
    return computeFreeEnergy(grid, rsu, coupling, true);
  }

  double computeFreeEnergy(const vector<double> &grid,
                           const vector<double> &rsu,
                           const double &coupling,
                           const bool normalize) {
    if (numUtil::largerThan(coupling, grid.back())) {
      parallelUtil::MPI::throwError(
          "The coupling parameter is out of range"
          " for the current grid, the free energy cannot be computed");
    }
    const Interpolator1D itp(grid, rsu);
    Integrator1D itg(1.0e-6);
    const FreeEnergy freeEnergy(coupling, itp, itg, normalize);
    return freeEnergy.get();
  }

  vector<double> computeRdf(const vector<double> &r,
                            const vector<double> &wvg,
                            const vector<double> &ssf) {
    assert(ssf.size() > 0 && wvg.size() > 0);
    const Interpolator1D itp(wvg, ssf);
    const int nr = r.size();
    vector<double> rdf(nr);
    Integrator1D itg(ItgType::DEFAULT, 1.0e-6);
    Integrator1D itgf(ItgType::FOURIER, 1.0e-6);
    for (int i = 0; i < nr; ++i) {
      const Rdf rdfTmp(r[i], wvg.back(), itp, itg, itgf);
      rdf[i] = rdfTmp.get();
    }
    return rdf;
  }

  // -----------------------------------------------------------------
  // InternalEnergy class
  // -----------------------------------------------------------------

  double InternalEnergy::ssf(const double &y) const { return ssfi.eval(y); }

  double InternalEnergy::integrand(const double &y) const { return ssf(y) - 1; }

  double InternalEnergy::get() const {
    auto func = [&](const double &y) -> double { return integrand(y); };
    itg.compute(func, ItgParam(yMin, yMax));
    return itg.getSolution() / (M_PI * rs * lambda);
  }

  // -----------------------------------------------------------------
  // FreeEnergy class
  // -----------------------------------------------------------------

  double FreeEnergy::get() const {
    auto func = [&](const double &y) -> double { return rsui.eval(y); };
    itg.compute(func, ItgParam(0.0, rs));
    if (normalize) {
      return (rs == 0.0) ? -numUtil::Inf : itg.getSolution() / rs / rs;
    };
    return itg.getSolution();
  }

  // -----------------------------------------------------------------
  // Rdf class
  // -----------------------------------------------------------------

  double Rdf::ssf(const double &y) const { return ssfi.eval(y); }

  double Rdf::integrand(const double &y) const {
    if (y > cutoff) return 0;
    const double yssf = y * (ssf(y) - 1);
    return (r == 0.0) ? y * yssf : yssf;
  }

  double Rdf::get() const {
    auto func = [&](const double &y) -> double { return integrand(y); };
    if (r == 0) {
      itg.compute(func, ItgParam(0.0, cutoff));
      return 1 + 1.5 * itg.getSolution();
    } else {
      itgf.compute(func, ItgParam(r));
      return 1 + 1.5 * itgf.getSolution() / r;
    }
  }

} // namespace thermoUtil

namespace parallelUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods for MPI distributed memory parallelism
  // -----------------------------------------------------------------

  namespace MPI {

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
      MPI_Allreduce(&myNumber,
                    &globalMininumNumber,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    MPICommunicator);
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

  } // namespace MPI

} // namespace parallelUtil
