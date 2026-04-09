#ifndef TESTS_VS_TEST_DOUBLES_HPP
#define TESTS_VS_TEST_DOUBLES_HPP

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "fixtures/input_builders.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"

namespace vsTestDoubles {

class FakeVSWorker : public VSWorker {
public:
  FakeVSWorker(double coupling,
               double degeneracy,
               double qadder,
               double uint_value,
               double chem_pot,
               double error_value,
               int tag)
      : coupling_(coupling),
        degeneracy_(degeneracy),
        qadder_(qadder),
        uint_(uint_value),
        chemical_potential_(chem_pot),
        error_(error_value),
        lfc_(3, 2),
        idr_(3, 2),
        wvg_{0.5, 1.0, 1.5},
        ssf_(wvg_.size(), 1.0),
        sdr_{static_cast<double>(tag)} {
    lfc_.fill(static_cast<double>(tag));
    idr_.fill(-static_cast<double>(tag));
  }

  const Vector2D &getLfc() const override { return lfc_; }
  const std::vector<double> &getWvg() const override { return wvg_; }
  const std::vector<double> &getSsf() const override { return ssf_; }
  const Vector2D &getIdr() const override { return idr_; }
  std::vector<double> getSdr() const override { return sdr_; }
  double getUInt() const override { return uint_; }
  double getChemicalPotential() const override { return chemical_potential_; }
  double getQAdder() const override { return qadder_; }
  double getCoupling() const override { return coupling_; }
  double getDegeneracy() const override { return degeneracy_; }

  void init() override { ++init_calls; }

  void initialGuess() override { ++initial_guess_calls; }

  void computeSsf() override {
    ++compute_ssf_calls;
    std::fill(ssf_.begin(), ssf_.end(), 1.0);
  }

  void computeLfc() override {
    ++compute_lfc_calls;
    for (size_t i = 0; i < lfc_.size(0); ++i) {
      for (size_t j = 0; j < lfc_.size(1); ++j) {
        lfc_(i, j) = static_cast<double>(i + j + 1);
      }
    }
  }

  void applyLfcDiff(const Vector2D &v) override {
    ++apply_lfc_diff_calls;
    last_lfc_diff = v;
    lfc_.diff(v);
  }

  double computeError() const override { return error_; }

  void updateSolution() override { ++update_solution_calls; }

  void setError(double value) { error_ = value; }

  int init_calls = 0;
  int initial_guess_calls = 0;
  int compute_ssf_calls = 0;
  int compute_lfc_calls = 0;
  int apply_lfc_diff_calls = 0;
  int update_solution_calls = 0;
  Vector2D last_lfc_diff;

private:
  double coupling_;
  double degeneracy_;
  double qadder_;
  double uint_;
  double chemical_potential_;
  double error_;
  Vector2D lfc_;
  Vector2D idr_;
  std::vector<double> wvg_;
  std::vector<double> ssf_;
  std::vector<double> sdr_;
};

class FakeVSManager : public VSManager {
public:
  explicit FakeVSManager(std::shared_ptr<const VSStlsInput> in)
      : in_(std::move(in)) {
    buildWorkersAndGridMetadata();
    setupDerivativeData();
  }

  int compute() override {
    init();
    initialGuess();
    computeSsf();
    computeLfc();
    updateSolution();
    return 0;
  }

  const FakeVSWorker &worker(GridPoint p) const {
    return *dynamic_cast<const FakeVSWorker *>(workers[p.toIndex()].get());
  }

  FakeVSWorker &worker(GridPoint p) {
    return *dynamic_cast<FakeVSWorker *>(workers[p.toIndex()].get());
  }

  void setCenterError(double value) { worker(GridPoints::CENTER).setError(value); }

protected:
  const VSInput &inVS() const override {
    return *static_cast<const VSInput *>(in_.get());
  }

  const Input &inScheme() const override {
    return *static_cast<const Input *>(in_.get());
  }

private:
  void buildWorkersAndGridMetadata() {
    const double drs = in_->getCouplingResolution();
    const double dtheta = in_->getDegeneracyResolution();
    const double rs0 = std::max(in_->getCoupling(), drs);
    const double theta0 = std::max(in_->getDegeneracy(), dtheta);

    for (const auto theta_off : {GridPoint::Theta::DOWN,
                                 GridPoint::Theta::CENTER,
                                 GridPoint::Theta::UP}) {
      const double theta = theta0 + static_cast<int>(theta_off) * dtheta;
      for (const auto rs_off :
           {GridPoint::Rs::DOWN, GridPoint::Rs::CENTER, GridPoint::Rs::UP}) {
        const double rs = rs0 + static_cast<int>(rs_off) * drs;
        const GridPoint p{rs_off, theta_off};
        const size_t idx = p.toIndex();
        rsValues[idx] = rs;
        thetaValues[idx] = theta;
        workers[idx] = std::make_unique<FakeVSWorker>(
            rs,
            theta,
            rs,
            10.0 + static_cast<double>(idx),
            20.0 + static_cast<double>(idx),
            0.01 * static_cast<double>(idx),
            static_cast<int>(idx));
      }
    }
  }

  std::shared_ptr<const VSStlsInput> in_;
};

class FakeVSBase : public VSBase {
public:
  explicit FakeVSBase(std::shared_ptr<const VSStlsInput> in)
      : in_(std::move(in)), manager_(in_) {
    setRsGrid();
    setFxcIntegrand();
  }

  const VSInput &in() const override {
    return *static_cast<const VSInput *>(in_.get());
  }

  const Input &inScheme() const override {
    return *static_cast<const Input *>(in_.get());
  }

  VSManager &grid() override { return manager_; }

  const VSManager &grid() const override { return manager_; }

  void resetGridStorage() {
    rsGrid.clear();
    fxcIntegrand.clear();
  }

  void callSetRsGrid() { setRsGrid(); }

  void callSetFxcIntegrand() { setFxcIntegrand(); }

  void callUpdateFxcIntegrand() { updateFxcIntegrand(); }

  int callRunGrid() { return runGrid(); }

  FakeVSManager &manager() { return manager_; }

  const FakeVSManager &manager() const { return manager_; }

private:
  std::shared_ptr<const VSStlsInput> in_;
  FakeVSManager manager_;
};

inline std::shared_ptr<VSStlsInput> makeVsInput(const double coupling = 1.0,
                                                const double degeneracy = 1.0) {
  auto in = std::make_shared<VSStlsInput>();
  testFixtures::configureIterationInput(
      *in, "VSSTLS", dimensionsUtil::Dimension::D3, coupling, degeneracy, 2);
  in->setCouplingResolution(0.5);
  in->setDegeneracyResolution(0.2);
  in->setAlphaGuess({0.2, 0.8});
  in->setErrMinAlpha(1.0e-6);
  in->setNIterAlpha(8);
  return in;
}

} // namespace vsTestDoubles

#endif
