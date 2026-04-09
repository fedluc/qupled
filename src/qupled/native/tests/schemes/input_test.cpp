#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "schemes/input.hpp"

TEST(InputTest, BaseInputSettersStoreValidValues) {
  Input in;
  in.setCoupling(2.0);
  in.setDegeneracy(0.5);
  in.setDimension(dimensionsUtil::Dimension::D2);
  in.setInt2DScheme("full");
  in.setIntError(1.0e-6);
  in.setNThreads(4);
  in.setTheory("STLS");
  in.setChemicalPotentialGuess({-1.0, 2.0});
  in.setNMatsubara(32);
  in.setWaveVectorGridRes(0.05);
  in.setWaveVectorGridCutoff(8.0);
  in.setFrequencyCutoff(20.0);

  EXPECT_DOUBLE_EQ(in.getCoupling(), 2.0);
  EXPECT_DOUBLE_EQ(in.getDegeneracy(), 0.5);
  EXPECT_EQ(in.getDimension(), dimensionsUtil::Dimension::D2);
  EXPECT_EQ(in.getInt2DScheme(), "full");
  EXPECT_DOUBLE_EQ(in.getIntError(), 1.0e-6);
  EXPECT_EQ(in.getNThreads(), 4);
  EXPECT_EQ(in.getTheory(), "STLS");
  EXPECT_TRUE(in.isClassic());
  EXPECT_EQ(in.getChemicalPotentialGuess(), std::vector<double>({-1.0, 2.0}));
  EXPECT_EQ(in.getNMatsubara(), 32);
  EXPECT_DOUBLE_EQ(in.getWaveVectorGridRes(), 0.05);
  EXPECT_DOUBLE_EQ(in.getWaveVectorGridCutoff(), 8.0);
  EXPECT_DOUBLE_EQ(in.getFrequencyCutoff(), 20.0);
}

TEST(InputTest, BaseInputSettersRejectInvalidValues) {
  Input in;
  EXPECT_THROW(in.setCoupling(-1.0), std::runtime_error);
  EXPECT_THROW(in.setDegeneracy(-1.0), std::runtime_error);
  EXPECT_THROW(in.setInt2DScheme("bad"), std::runtime_error);
  EXPECT_THROW(in.setIntError(0.0), std::runtime_error);
  EXPECT_THROW(in.setNThreads(0), std::runtime_error);
  EXPECT_THROW(in.setTheory("unknown"), std::runtime_error);
  EXPECT_THROW(in.setChemicalPotentialGuess({1.0, 1.0}), std::runtime_error);
  EXPECT_THROW(in.setChemicalPotentialGuess({1.0}), std::runtime_error);
  EXPECT_THROW(in.setNMatsubara(-1), std::runtime_error);
  EXPECT_THROW(in.setWaveVectorGridRes(0.0), std::runtime_error);
  EXPECT_THROW(in.setWaveVectorGridCutoff(0.0), std::runtime_error);
  EXPECT_THROW(in.setFrequencyCutoff(0.0), std::runtime_error);
}

TEST(InputTest, IterationInputValidationAndGuessConsistencyWork) {
  IterationInput in;
  in.setErrMin(1.0e-5);
  in.setMixingParameter(0.5);
  in.setNIter(20);
  const Guess guess{
      .wvg = {1.0, 2.0},
      .ssf = {0.9, 1.1},
      .lfc = Vector2D(2, 1),
  };
  in.setGuess(guess);

  EXPECT_DOUBLE_EQ(in.getErrMin(), 1.0e-5);
  EXPECT_DOUBLE_EQ(in.getMixingParameter(), 0.5);
  EXPECT_EQ(in.getNIter(), 20);
  EXPECT_EQ(in.getGuess().wvg.size(), 2u);

  EXPECT_THROW(in.setErrMin(0.0), std::runtime_error);
  EXPECT_THROW(in.setMixingParameter(-0.1), std::runtime_error);
  EXPECT_THROW(in.setMixingParameter(1.1), std::runtime_error);
  EXPECT_THROW(in.setNIter(-1), std::runtime_error);
  EXPECT_THROW(in.setGuess(Guess{{1.0}, {1.0, 2.0}, Vector2D(1, 1)}),
               std::runtime_error);
}

TEST(InputTest, IetAndVsInputsValidateMappingsAndGridShapes) {
  IetInput iet;
  iet.setMapping("linear");
  EXPECT_EQ(iet.getMapping(), "linear");
  EXPECT_THROW(iet.setMapping("invalid"), std::runtime_error);

  VSInput vs;
  vs.setAlphaGuess({0.1, 1.0});
  vs.setCouplingResolution(0.01);
  vs.setDegeneracyResolution(0.02);
  vs.setErrMinAlpha(1.0e-6);
  vs.setNIterAlpha(30);

  VSInput::FreeEnergyIntegrand valid_integrand;
  valid_integrand.grid = {0.2, 0.4, 0.6};
  valid_integrand.integrand = {{1.0, 1.5, 2.0}, {2.0, 2.5, 3.0}};
  vs.setFreeEnergyIntegrand(valid_integrand);
  EXPECT_EQ(vs.getFreeEnergyIntegrand().grid.size(), 3u);

  EXPECT_THROW(vs.setAlphaGuess({1.0, 1.0}), std::runtime_error);
  EXPECT_THROW(vs.setCouplingResolution(0.0), std::runtime_error);
  EXPECT_THROW(vs.setDegeneracyResolution(0.0), std::runtime_error);
  EXPECT_THROW(vs.setErrMinAlpha(0.0), std::runtime_error);
  EXPECT_THROW(vs.setNIterAlpha(-1), std::runtime_error);

  VSInput::FreeEnergyIntegrand bad_integrand;
  bad_integrand.grid = {0.2, 0.4};
  bad_integrand.integrand = {{1.0}, {2.0, 3.0}};
  EXPECT_THROW(vs.setFreeEnergyIntegrand(bad_integrand), std::runtime_error);
}
