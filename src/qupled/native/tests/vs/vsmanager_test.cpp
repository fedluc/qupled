#include <gtest/gtest.h>

#include "vs/grid_point.hpp"
#include "vs/test_doubles.hpp"

using namespace vsTestDoubles;

TEST(VSManagerTest, SetAlphaStoresProvidedValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  manager.setAlpha(0.35);

  EXPECT_DOUBLE_EQ(manager.getAlpha(), 0.35);
}

TEST(VSManagerTest, ComputeInitializesEachWorkerExactlyOnce) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();
  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).init_calls, 1);
}

TEST(VSManagerTest, ComputeRunsInitialGuessOnWorkers) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).initial_guess_calls, 1);
}

TEST(VSManagerTest, ComputeRunsSsfStepOnWorkers) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).compute_ssf_calls, 1);
}

TEST(VSManagerTest, ComputeRunsLfcStepOnWorkers) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).compute_lfc_calls, 1);
}

TEST(VSManagerTest, ComputeAppliesLfcCorrectionOnWorkers) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).apply_lfc_diff_calls, 1);
}

TEST(VSManagerTest, ComputeRunsSolutionUpdateOnWorkers) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setAlpha(0.5);

  manager.compute();

  EXPECT_EQ(manager.worker(GridPoints::CENTER).update_solution_calls, 1);
}

TEST(VSManagerTest, GetErrorUsesCenterWorkerError) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  manager.setCenterError(0.123);

  EXPECT_DOUBLE_EQ(manager.getError(), 0.123);
}

TEST(VSManagerTest, GetSsfReturnsSelectedGridPointData) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  const auto &ssf = manager.getSsf(GridPoints::RS_UP_THETA_UP);

  ASSERT_EQ(ssf.size(), 3u);
  EXPECT_DOUBLE_EQ(ssf[0], 1.0);
}

TEST(VSManagerTest, GetLfcReturnsSelectedGridPointData) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  const auto &lfc = manager.getLfc(GridPoints::CENTER);

  EXPECT_EQ(lfc.size(0), 3u);
  EXPECT_EQ(lfc.size(1), 2u);
  EXPECT_DOUBLE_EQ(lfc(0, 0), 4.0);
}

TEST(VSManagerTest, GetWvgReturnsSelectedGridPointData) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  const auto &wvg = manager.getWvg(GridPoints::CENTER);

  ASSERT_EQ(wvg.size(), 3u);
  EXPECT_DOUBLE_EQ(wvg[0], 0.5);
}

TEST(VSManagerTest, GetIdrReturnsSelectedGridPointData) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  const auto &idr = manager.getIdr(GridPoints::RS_UP_THETA);

  EXPECT_EQ(idr.size(0), 3u);
  EXPECT_EQ(idr.size(1), 2u);
  EXPECT_DOUBLE_EQ(idr(0, 0), -5.0);
}

TEST(VSManagerTest, GetSdrReturnsSelectedGridPointData) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  const auto sdr = manager.getSdr(GridPoints::RS_DOWN_THETA_UP);

  ASSERT_EQ(sdr.size(), 1u);
  EXPECT_DOUBLE_EQ(sdr[0], 6.0);
}

TEST(VSManagerTest, GetCouplingReturnsSelectedGridPointValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_DOUBLE_EQ(manager.getCoupling(GridPoints::CENTER), 1.0);
}

TEST(VSManagerTest, GetDegeneracyReturnsSelectedGridPointValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_DOUBLE_EQ(manager.getDegeneracy(GridPoints::CENTER), 1.0);
}

TEST(VSManagerTest, GetUIntReturnsSelectedGridPointValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_DOUBLE_EQ(manager.getUInt(GridPoints::RS_UP_THETA), 15.0);
}

TEST(VSManagerTest, GetChemicalPotentialReturnsSelectedGridPointValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_DOUBLE_EQ(manager.getChemicalPotential(GridPoints::RS_UP_THETA), 25.0);
}

TEST(VSManagerTest, GetQAdderReturnsSelectedGridPointValue) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_DOUBLE_EQ(manager.getQAdder(GridPoints::CENTER), 1.0);
}

TEST(VSManagerTest, GetFxcIntegrandValueReturnsZeroForUnitStructureFactor) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);

  EXPECT_NEAR(manager.getFxcIntegrandValue(GridPoints::CENTER), 0.0, 1.0e-12);
}

#ifndef NDEBUG
TEST(VSManagerTest, ComputeDiesWhenAlphaWasNotSet) {
  const auto in = makeVsInput();
  FakeVSManager manager(in);
  EXPECT_DEATH((void)manager.compute(), "");
}
#endif
