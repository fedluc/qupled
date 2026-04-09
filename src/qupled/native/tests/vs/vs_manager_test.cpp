#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <vector>

#include "schemes/input.hpp"
#include "vs/vsmanager.hpp"

namespace {

class FakeWorker : public VSWorker {
public:
  FakeWorker()
      : lfc_(3, 1),
        idr_(3, 1),
        wvg_{0.5, 1.0, 1.5},
        ssf_{1.0, 1.0, 1.0} {}

  const Vector2D &getLfc() const override { return lfc_; }
  const std::vector<double> &getWvg() const override { return wvg_; }
  const std::vector<double> &getSsf() const override { return ssf_; }
  const Vector2D &getIdr() const override { return idr_; }
  std::vector<double> getSdr() const override { return {0.0, 0.0, 0.0}; }
  double getUInt() const override { return 0.0; }
  double getChemicalPotential() const override { return 0.0; }
  double getQAdder() const override { return 0.0; }
  double getCoupling() const override { return 1.0; }
  double getDegeneracy() const override { return 1.0; }

  void init() override {}
  void initialGuess() override {}
  void computeSsf() override {}
  void computeLfc() override {
    ++computeLfcCalls;
    for (size_t i = 0; i < wvg_.size(); ++i) {
      lfc_(i, 0) = static_cast<double>(i + 1);
    }
  }
  void applyLfcDiff(const Vector2D &v) override {
    ++applyDiffCalls;
    lfc_.diff(v);
  }
  double computeError() const override { return 0.0; }
  void updateSolution() override {}

  int computeLfcCalls = 0;
  int applyDiffCalls = 0;

private:
  Vector2D lfc_;
  Vector2D idr_;
  std::vector<double> wvg_;
  std::vector<double> ssf_;
};

class TestVSManager : public VSManager {
public:
  TestVSManager() {
    in_.setTheory("VSSTLS");
    in_.setDimension(dimensionsUtil::Dimension::D3);
    in_.setCoupling(1.0);
    in_.setDegeneracy(1.0);
    in_.setWaveVectorGridRes(0.5);
    in_.setCouplingResolution(0.2);
    in_.setDegeneracyResolution(0.2);
    for (size_t i = 0; i < N; ++i) {
      workers[i] = std::make_unique<FakeWorker>();
      rsValues[i] = 1.0;
      thetaValues[i] = 1.0;
    }
    setupDerivativeData();
  }

  int compute() override { return 0; }

  void setAlphaForTest(double a) { alpha = a; }
  void computeLfcForTest() { computeLfc(); }

  bool isRsCentered(GridPoint p) const {
    return rsDerivData[p.toIndex()].type == DerivativeData::Type::CENTERED;
  }

  bool isRsForward(GridPoint p) const {
    return rsDerivData[p.toIndex()].type == DerivativeData::Type::FORWARD;
  }

  bool isThetaBackward(GridPoint p) const {
    return thetaDerivData[p.toIndex()].type == DerivativeData::Type::BACKWARD;
  }

  const FakeWorker &worker(GridPoint p) const {
    return *static_cast<FakeWorker *>(workers[p.toIndex()].get());
  }

private:
  VSStlsInput in_;

  const VSInput &inVS() const override { return in_; }
  const Input &inScheme() const override { return in_; }
};

} // namespace

TEST(VSManagerTest, DerivativeStencilMetadataUsesExpectedNeighbors) {
  const TestVSManager mgr;
  EXPECT_TRUE(mgr.isRsCentered(GridPoints::CENTER));
  EXPECT_TRUE(mgr.isRsForward(GridPoints::RS_DOWN_THETA));
  EXPECT_TRUE(mgr.isThetaBackward(GridPoints::RS_THETA_UP));
}

TEST(VSManagerTest, ComputeLfcRunsWorkersAndAppliesDerivativeCorrection) {
  TestVSManager mgr;
  mgr.setAlphaForTest(0.5);
  mgr.computeLfcForTest();

  for (const auto p : {GridPoints::RS_DOWN_THETA_DOWN,
                       GridPoints::RS_THETA_DOWN,
                       GridPoints::RS_UP_THETA_DOWN,
                       GridPoints::RS_DOWN_THETA,
                       GridPoints::CENTER,
                       GridPoints::RS_UP_THETA,
                       GridPoints::RS_DOWN_THETA_UP,
                       GridPoints::RS_THETA_UP,
                       GridPoints::RS_UP_THETA_UP}) {
    const auto &w = mgr.worker(p);
    EXPECT_EQ(w.computeLfcCalls, 1);
    EXPECT_EQ(w.applyDiffCalls, 1);
  }
}
