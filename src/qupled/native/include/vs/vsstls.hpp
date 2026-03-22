#ifndef VS_VSSTLS_HPP
#define VS_VSSTLS_HPP

#include "input.hpp"
#include "stls.hpp"
#include "vs/state_point_grid.hpp"
#include "vs/vsbase.hpp"

// -----------------------------------------------------------------
// VSStlsWorker: implements IVSWorker for STLS-based VS workers
// -----------------------------------------------------------------

class VSStlsWorker : public IVSWorker, public Stls {
public:
  using InputType = VSStlsInput;

  VSStlsWorker(const std::shared_ptr<const VSStlsInput> &in,
               bool verbose,
               GridPoint /* unused */)
      : Stls(in, verbose) {}

  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  void computeBaseLfc()                override { Stls::computeLfc(); }

  const Vector2D &            getLfc() const override { return lfc; }
  const std::vector<double> & getWvg() const override { return wvg; }
  const std::vector<double> & getSsf() const override { return ssf; }

  void   workerInit()              override { HF::init(); }
  void   workerInitialGuess()      override { Stls::initialGuess(); }
  void   workerComputeSsf()        override { Stls::computeSsf(); }
  double workerComputeError() const override { return Stls::computeError(); }
  void   workerUpdateSolution()    override { Stls::updateSolution(); }
};

// -----------------------------------------------------------------
// StatePointGridVSStls: manages 9 VSStlsWorkers, inherits Stls loop
// -----------------------------------------------------------------

class StatePointGridVSStls : public StatePointGridBase, public Stls {
public:
  explicit StatePointGridVSStls(const std::shared_ptr<const VSStlsInput> &in);

protected:
  void   init()              override;
  void   computeLfc()        override;
  void   computeSsf()        override;
  double computeError() const override;
  void   updateSolution()    override;
  void   initialGuess()      override;
};

// -----------------------------------------------------------------
// VSStls: VS scheme based on Stls structural properties
// -----------------------------------------------------------------

class VSStls : public VSBase {
public:
  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in);
  using VSBase::compute;

  const std::vector<double> &getSsf() const;
  const Vector2D &            getLfc() const;
  const std::vector<double> &getWvg() const;
  const Vector2D &            getIdr() const;
  std::vector<double>         getSdr() const;
  double                      getUInt() const;
  double                      getError() const;

private:
  std::shared_ptr<const VSStlsInput> inPtr;
  StatePointGridVSStls               grid;

  const VSInput &in()       const override;
  const Input &  inScheme() const override;

  int    runGrid()                               override;
  double getCoupling(GridPoint p)          const override;
  double getDegeneracy(GridPoint p)        const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  double computeQRaw(GridPoint p)          const override;
};

#endif
