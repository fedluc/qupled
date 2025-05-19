#ifndef IET_HPP
#define IET_HPP

#include "input.hpp"
#include "logger.hpp"
#include "numerics.hpp"
#include "stls.hpp"
#include "vector2D.hpp"

// -----------------------------------------------------------------
// Solver for the STLS scheme
// -----------------------------------------------------------------

class Iet : public Logger {

public:

  // Constructors
  explicit Iet(const std::shared_ptr<const IetInput> &in_,
               const std::vector<double> &wvg_)
      : Logger(true),
        inPtr(in_),
        wvg(wvg_),
        bf(wvg_.size()) {}
  // Initialize the iet object
  void init();
  // Getters
  const std::vector<double> &getBf() const { return bf; }
  // Read initial guess from input
  bool initialGuessFromInput(Vector2D &lfc);

private:

  // Input parameters
  const std::shared_ptr<const IetInput> inPtr;
  // Access input pointer
  const IetInput &in() const { return *inPtr; }
  // Cast the input member to an IterationInput type
  const IterationInput &inRpa() const {
    return *StlsUtil::dynamic_pointer_cast<IetInput, IterationInput>(inPtr);
  }
  // Wave vector grid
  const std::vector<double> wvg;
  // Bridge function
  std::vector<double> bf;
  // Compute bridge function
  void computeBf();
};

namespace IetUtil {

  class BridgeFunction {

  public:

    // Constructor
    BridgeFunction(const std::string &theory_,
                   const std::string &mapping_,
                   const double &rs_,
                   const double &Theta_,
                   const double &x_,
                   std::shared_ptr<Integrator1D> itg_)
        : theory(theory_),
          mapping(mapping_),
          rs(rs_),
          Theta(Theta_),
          x(x_),
          itg(itg_) {}
    // Get result of the integration
    double get() const;

  private:

    // Theory to be solved
    const std::string theory;
    // Iet mapping
    const std::string mapping;
    // Coupling parameter
    const double rs;
    // Degeneracy parameter
    const double Theta;
    // Wave vector
    const double x;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Constant for unit conversion
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
    // Hypernetted-chain bridge function
    double hnc() const;
    // Ichimaru bridge function
    double ioi() const;
    // Lucco Castello and Tolias bridge function
    double lct() const;
    double lctIntegrand(const double &r, const double &Gamma) const;
    // Coupling parameter to compute the bridge function
    double couplingParameter() const;
  };

} // namespace IetUtil

#endif
