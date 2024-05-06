#ifndef UTIL_HPP
#define UTIL_HPP

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

class Interpolator1D;
class Integrator1D;

// -----------------------------------------------------------------
// Utility functions to compute thermodynamic properties
// -----------------------------------------------------------------

namespace thermoUtil {

  double computeInternalEnergy(const std::vector<double> &wvg,
                               const std::vector<double> &ssf,
                               const double &coupling);

  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling);

  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling,
                           const bool normalize);

  std::vector<double> computeRdf(const std::vector<double> &r,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &ssf);

  // --- Class for internal energy calculation ---
  class InternalEnergy {

  public:

    // Constructor
    InternalEnergy(const double &rs_,
                   const double &yMin_,
                   const double &yMax_,
                   const Interpolator1D &ssfi_,
                   Integrator1D &itg_)
        : rs(rs_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_),
          ssfi(ssfi_) {}
    // Get result of integration
    double get() const;

  private:

    // Coupling parameter
    const double rs;
    // Integration limits
    const double yMin;
    const double yMax;
    // Integrator object
    Integrator1D &itg;
    // Static structure factor interpolator
    const Interpolator1D &ssfi;
    // Integrand
    double integrand(const double &y) const;
    // Compute static structure factor
    double ssf(const double &y) const;
    // Constant for unit conversion
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
  };

  // --- Class for free energy calculation ---
  class FreeEnergy {

  public:

    // Constructor
    FreeEnergy(const double &rs_,
               const Interpolator1D &rsui_,
               Integrator1D &itg_,
               const bool normalize_)
        : rs(rs_),
          itg(itg_),
          rsui(rsui_),
          normalize(normalize_) {}
    // Get result of integration
    double get() const;

  private:

    // Coupling parameter
    const double rs;
    // Integrator object
    Integrator1D &itg;
    // Integrand interpolator (the integrand is given by rs * u)
    const Interpolator1D &rsui;
    // Integrand
    double integrand(const double y) const;
    // Flag marking whether the free energy should be normalized with rs^2
    const bool normalize;
  };

  // --- Class for radial distribution function calculation ---
  class Rdf {

  public:

    // Constructor
    Rdf(const double &r_,
        const double &cutoff_,
        const Interpolator1D &ssfi_,
        Integrator1D &itg_,
        Integrator1D &itgf_)
        : r(r_),
          cutoff(cutoff_),
          itgf(itgf_),
          itg(itg_),
          ssfi(ssfi_) {}
    // Get result of integration
    double get() const;

  private:

    // Spatial position
    const double r;
    // Cutoff in the wave-vector grid
    const double cutoff;
    // Fourier Integrator object
    Integrator1D &itgf;
    // Integrator object
    Integrator1D &itg;
    // Static structure factor interpolator
    const Interpolator1D &ssfi;
    // Integrand
    double integrand(const double &y) const;
    // Compute static structure factor
    double ssf(const double &y) const;
  };

} // namespace thermoUtil

// -----------------------------------------------------------------
// Utility functions to manipulate binary files
// -----------------------------------------------------------------

namespace binUtil {

  template <typename T> void writeNum(std::ofstream &file, const T &num) {
    file.write(reinterpret_cast<const char *>(&num), sizeof(num));
  }

  template <typename T>
  void writeDataToBinary(std::ofstream &file, const double &data) {
    writeNum<double>(file, data);
  }

  template <typename T>
  void writeDataToBinary(std::ofstream &file, const int &data) {
    writeNum<int>(file, data);
  }

  template <class T>
  void writeDataToBinary(std::ofstream &file, const T &data) {
    for (auto &el : data) {
      writeDataToBinary<decltype(el)>(file, el);
    }
  }

  template <typename T> void readNum(std::ifstream &file, T &num) {
    file.read((char *)&num, sizeof(T));
  }

  template <typename T>
  void readDataFromBinary(std::ifstream &file, double &data) {
    readNum<double>(file, data);
  }

  template <typename T>
  void readDataFromBinary(std::ifstream &file, int &data) {
    readNum<int>(file, data);
  }

  template <class T> void readDataFromBinary(std::ifstream &file, T &data) {
    for (auto &el : data) {
      readDataFromBinary<decltype(el)>(file, el);
    }
  }

} // namespace binUtil

// -------------------------------------------------------------------
// Utility functions to handle parallel calculations with OMP and MPI
// -------------------------------------------------------------------

namespace parallelUtil {

  // --- MPI for distributed memory parallelism ---
  namespace MPI {

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

  } // namespace MPI

} // namespace parallelUtil

#endif
