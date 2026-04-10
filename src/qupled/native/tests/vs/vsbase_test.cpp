#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "vs/grid_point.hpp"
#include "vs/test_doubles.hpp"

using namespace vsTestDoubles;

TEST(VSBaseTest, SetRsGridBuildsExpectedGridValues) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &grid = base.getFreeEnergyGrid();

  ASSERT_EQ(grid.size(), 4u);
  EXPECT_DOUBLE_EQ(grid[0], 0.0);
  EXPECT_DOUBLE_EQ(grid[1], 0.5);
  EXPECT_DOUBLE_EQ(grid[2], 1.0);
  EXPECT_DOUBLE_EQ(grid[3], 1.5);
}

TEST(VSBaseTest, SetRsGridRejectsCouplingNotMultipleOfResolution) {
  const auto in = makeVsInput(1.1, 1.0);
  EXPECT_THROW((FakeVSBase(in)), std::runtime_error);
}

TEST(VSBaseTest, SetFxcIntegrandInitializesThreeThetaRows) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &integrand = base.getFreeEnergyIntegrand();

  ASSERT_EQ(integrand.size(), 3u);
  EXPECT_EQ(integrand[0].size(), base.getFreeEnergyGrid().size());
  EXPECT_EQ(integrand[1].size(), base.getFreeEnergyGrid().size());
  EXPECT_EQ(integrand[2].size(), base.getFreeEnergyGrid().size());
}

TEST(VSBaseTest, SetFxcIntegrandUsesInfinityAsDefaultSentinel) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_TRUE(std::isinf(integrand[0][0]));
  EXPECT_TRUE(std::isinf(integrand[1][0]));
  EXPECT_TRUE(std::isinf(integrand[2][0]));
}

TEST(VSBaseTest, SetFxcIntegrandLoadsPrecomputedDataWhenProvided) {
  const auto in = makeVsInput(1.0, 1.0);
  VSInput::FreeEnergyIntegrand preload;
  preload.grid = {0.0, 0.5, 1.0, 1.5};
  preload.integrand = {{10.0, 11.0, 12.0, 13.0},
                       {20.0, 21.0, 22.0, 23.0},
                       {30.0, 31.0, 32.0, 33.0}};
  in->setFreeEnergyIntegrand(preload);

  FakeVSBase base(in);
  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_DOUBLE_EQ(integrand[0][2], 12.0);
  EXPECT_DOUBLE_EQ(integrand[1][2], 22.0);
  EXPECT_DOUBLE_EQ(integrand[2][2], 32.0);
}

TEST(VSBaseTest, UpdateFxcIntegrandUpdatesDownThetaColumnNeighbors) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  base.callUpdateFxcIntegrand();
  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_NEAR(integrand[0][1], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[0][2], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[0][3], 0.0, 1.0e-12);
}

TEST(VSBaseTest, UpdateFxcIntegrandUpdatesCenterThetaColumnNeighbors) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  base.callUpdateFxcIntegrand();
  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_NEAR(integrand[1][1], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[1][2], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[1][3], 0.0, 1.0e-12);
}

TEST(VSBaseTest, UpdateFxcIntegrandUpdatesUpThetaColumnNeighbors) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  base.callUpdateFxcIntegrand();
  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_NEAR(integrand[2][1], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[2][2], 0.0, 1.0e-12);
  EXPECT_NEAR(integrand[2][3], 0.0, 1.0e-12);
}

TEST(VSBaseTest, UpdateFxcIntegrandLeavesUnrelatedGridEntriesUntouched) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  base.callUpdateFxcIntegrand();
  const auto &integrand = base.getFreeEnergyIntegrand();

  EXPECT_TRUE(std::isinf(integrand[0][0]));
  EXPECT_TRUE(std::isinf(integrand[1][0]));
  EXPECT_TRUE(std::isinf(integrand[2][0]));
}

TEST(VSBaseTest, RunGridReturnsSuccessStatus) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);
  base.manager().setAlpha(0.5);

  EXPECT_EQ(base.callRunGrid(), 0);
}

TEST(VSBaseTest, GetErrorReadsFromManager) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);
  base.manager().setCenterError(0.222);

  EXPECT_DOUBLE_EQ(base.getError(), 0.222);
}

TEST(VSBaseTest, GetSsfReturnsCenterGridPointDataForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &ssf = base.getSsf();

  ASSERT_EQ(ssf.size(), 3u);
  EXPECT_DOUBLE_EQ(ssf[0], 1.0);
}

TEST(VSBaseTest, GetLfcReturnsCenterGridPointDataForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &lfc = base.getLfc();

  EXPECT_EQ(lfc.size(0), 3u);
  EXPECT_EQ(lfc.size(1), 2u);
  EXPECT_DOUBLE_EQ(lfc(0, 0), 4.0);
}

TEST(VSBaseTest, GetWvgReturnsCenterGridPointDataForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &wvg = base.getWvg();

  ASSERT_EQ(wvg.size(), 3u);
  EXPECT_DOUBLE_EQ(wvg[0], 0.5);
}

TEST(VSBaseTest, GetIdrReturnsCenterGridPointDataForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto &idr = base.getIdr();

  EXPECT_EQ(idr.size(0), 3u);
  EXPECT_EQ(idr.size(1), 2u);
  EXPECT_DOUBLE_EQ(idr(0, 0), -4.0);
}

TEST(VSBaseTest, GetSdrReturnsCenterGridPointDataForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);

  const auto sdr = base.getSdr();

  ASSERT_EQ(sdr.size(), 1u);
  EXPECT_DOUBLE_EQ(sdr[0], 4.0);
}

TEST(VSBaseTest, GetUIntUsesCenterOutputPointForFiniteState) {
  const auto in = makeVsInput(1.0, 1.0);
  FakeVSBase base(in);
  EXPECT_DOUBLE_EQ(base.getUInt(), 14.0);
}

TEST(VSBaseTest, GetUIntUsesRsDownThetaOutputPointAtZeroCoupling) {
  const auto in = makeVsInput(0.0, 1.0);
  FakeVSBase base(in);
  EXPECT_DOUBLE_EQ(base.getUInt(), 13.0);
}

TEST(VSBaseTest, GetUIntUsesRsThetaDownOutputPointAtZeroDegeneracy) {
  const auto in = makeVsInput(1.0, 0.0);
  FakeVSBase base(in);
  EXPECT_DOUBLE_EQ(base.getUInt(), 11.0);
}

TEST(VSBaseTest, GetUIntUsesRsDownThetaDownOutputPointAtDoubleZeroState) {
  const auto in = makeVsInput(0.0, 0.0);
  FakeVSBase base(in);
  EXPECT_DOUBLE_EQ(base.getUInt(), 10.0);
}

TEST(VSBaseTest, GetChemicalPotentialUsesMappedOutputPoint) {
  const auto in = makeVsInput(0.0, 0.0);
  FakeVSBase base(in);
  EXPECT_DOUBLE_EQ(base.getChemicalPotential(), 20.0);
}

TEST(VSBaseTest, ComputeConvergesAndReturnsSuccess) {
  const auto in = makeVsInput(1.0, 1.0);
  VSInput::FreeEnergyIntegrand preload;
  preload.grid = {0.0, 0.5, 1.0, 1.5};
  preload.integrand = {
      {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  in->setFreeEnergyIntegrand(preload);
  FakeVSBase base(in);

  EXPECT_EQ(base.compute(), 0);
  EXPECT_NEAR(base.getAlpha(), 1.0, 1.0e-8);
}
