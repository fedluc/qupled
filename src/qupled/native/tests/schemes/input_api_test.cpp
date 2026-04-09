#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "schemes/input.hpp"

TEST(InputApiTest, BaseSettersRoundTripAndTheoryFlags) {
  Input in;
  databaseUtil::DatabaseInfo db;
  db.name = "tmp.sqlite";
  db.blobStorage = "/tmp";
  db.runId = 7;
  db.runTableName = "runs";

  in.setCoupling(1.5);
  in.setDatabaseInfo(db);
  in.setDimension(dimensionsUtil::Dimension::D2);
  in.setDegeneracy(0.4);
  in.setInt2DScheme("segregated");
  in.setIntError(1.0e-7);
  in.setNThreads(2);
  in.setTheory("QSTLS");
  in.setChemicalPotentialGuess({-3.0, 2.0});
  in.setNMatsubara(6);
  in.setWaveVectorGridRes(0.1);
  in.setWaveVectorGridCutoff(3.0);
  in.setFrequencyCutoff(12.0);

  EXPECT_DOUBLE_EQ(in.getCoupling(), 1.5);
  EXPECT_EQ(in.getDatabaseInfo().runId, 7);
  EXPECT_EQ(in.getDimension(), dimensionsUtil::Dimension::D2);
  EXPECT_DOUBLE_EQ(in.getDegeneracy(), 0.4);
  EXPECT_EQ(in.getInt2DScheme(), "segregated");
  EXPECT_DOUBLE_EQ(in.getIntError(), 1.0e-7);
  EXPECT_EQ(in.getNThreads(), 2);
  EXPECT_EQ(in.getTheory(), "QSTLS");
  EXPECT_FALSE(in.isClassic());
  EXPECT_EQ(in.getChemicalPotentialGuess(), std::vector<double>({-3.0, 2.0}));
  EXPECT_EQ(in.getNMatsubara(), 6);
  EXPECT_DOUBLE_EQ(in.getWaveVectorGridRes(), 0.1);
  EXPECT_DOUBLE_EQ(in.getWaveVectorGridCutoff(), 3.0);
  EXPECT_DOUBLE_EQ(in.getFrequencyCutoff(), 12.0);
}

TEST(InputApiTest, SettersRejectInvalidValues) {
  Input in;
  EXPECT_THROW(in.setCoupling(-1.0), std::runtime_error);
  EXPECT_THROW(in.setDegeneracy(-0.1), std::runtime_error);
  EXPECT_THROW(in.setInt2DScheme("invalid"), std::runtime_error);
  EXPECT_THROW(in.setIntError(0.0), std::runtime_error);
  EXPECT_THROW(in.setNThreads(0), std::runtime_error);
  EXPECT_THROW(in.setTheory("NOT_A_THEORY"), std::runtime_error);
  EXPECT_THROW(in.setChemicalPotentialGuess({1.0}), std::runtime_error);
  EXPECT_THROW(in.setChemicalPotentialGuess({1.0, 1.0}), std::runtime_error);
  EXPECT_THROW(in.setNMatsubara(-1), std::runtime_error);
  EXPECT_THROW(in.setWaveVectorGridRes(0.0), std::runtime_error);
  EXPECT_THROW(in.setWaveVectorGridCutoff(0.0), std::runtime_error);
  EXPECT_THROW(in.setFrequencyCutoff(0.0), std::runtime_error);
}

TEST(InputApiTest, IterationIetQuantumAndVsInputsValidateContracts) {
  IterationInput iter;
  iter.setErrMin(1.0e-4);
  iter.setMixingParameter(0.5);
  iter.setNIter(3);
  Guess guess{.wvg = {0.0, 1.0}, .ssf = {0.2, 0.9}, .lfc = Vector2D(2, 1)};
  iter.setGuess(guess);
  EXPECT_DOUBLE_EQ(iter.getErrMin(), 1.0e-4);
  EXPECT_DOUBLE_EQ(iter.getMixingParameter(), 0.5);
  EXPECT_EQ(iter.getNIter(), 3);
  EXPECT_EQ(iter.getGuess().wvg.size(), 2u);
  EXPECT_THROW(iter.setErrMin(0.0), std::runtime_error);
  EXPECT_THROW(iter.setMixingParameter(1.1), std::runtime_error);
  EXPECT_THROW(iter.setNIter(-1), std::runtime_error);
  EXPECT_THROW(iter.setGuess(Guess{{1.0}, {1.0, 2.0}, Vector2D(1, 1)}),
               std::runtime_error);

  QstlsInput qin;
  qin.setFixedRunId(42);
  EXPECT_EQ(qin.getFixedRunId(), 42);

  IetInput iet;
  iet.setMapping("sqrt");
  EXPECT_EQ(iet.getMapping(), "sqrt");
  EXPECT_THROW(iet.setMapping("bad"), std::runtime_error);

  VSInput vs;
  vs.setAlphaGuess({0.2, 0.8});
  vs.setCouplingResolution(0.1);
  vs.setDegeneracyResolution(0.2);
  vs.setErrMinAlpha(1.0e-6);
  vs.setNIterAlpha(4);
  VSInput::FreeEnergyIntegrand f;
  f.grid = {0.4, 0.8};
  f.integrand = {{1.0, 2.0}, {3.0, 4.0}};
  vs.setFreeEnergyIntegrand(f);
  EXPECT_EQ(vs.getAlphaGuess(), std::vector<double>({0.2, 0.8}));
  EXPECT_EQ(vs.getNIterAlpha(), 4);
  EXPECT_EQ(vs.getFreeEnergyIntegrand().grid.size(), 2u);
  EXPECT_THROW(vs.setAlphaGuess({1.0, 1.0}), std::runtime_error);
  EXPECT_THROW(vs.setCouplingResolution(0.0), std::runtime_error);
  EXPECT_THROW(vs.setDegeneracyResolution(0.0), std::runtime_error);
  EXPECT_THROW(vs.setErrMinAlpha(0.0), std::runtime_error);
  EXPECT_THROW(vs.setNIterAlpha(-1), std::runtime_error);
}
