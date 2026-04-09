#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "fixtures/input_builders.hpp"
#include "schemes/vsstls.hpp"
#include "vs/grid_point.hpp"

TEST(VsStlsApiTest, ConstructorGuardRejectsUnsupported2D) {
  auto in = std::make_shared<VSStlsInput>();
  testFixtures::configureIterationInput(
      *in, "VSSTLS", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);
  in->setCouplingResolution(0.2);
  in->setDegeneracyResolution(0.2);
  in->setAlphaGuess({0.2, 0.8});
  in->setErrMinAlpha(1.0e-4);
  in->setNIterAlpha(3);

  EXPECT_THROW((VSStls(in)), std::runtime_error);
}

TEST(VsStlsApiTest, ConstructorAcceptsValid3DConfiguration) {
  auto in = std::make_shared<VSStlsInput>();
  testFixtures::configureIterationInput(
      *in, "VSSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  in->setCouplingResolution(0.2);
  in->setDegeneracyResolution(0.2);
  in->setAlphaGuess({0.2, 0.8});
  in->setErrMinAlpha(1.0e-4);
  in->setNIterAlpha(3);

  EXPECT_NO_THROW((VSStls(in)));
}

TEST(VsStlsApiTest, ManagerAndWorkerApiExposeGridContracts) {
  auto in = std::make_shared<VSStlsInput>();
  testFixtures::configureIterationInput(
      *in, "VSSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.7, 2);
  in->setCouplingResolution(0.2);
  in->setDegeneracyResolution(0.2);
  in->setAlphaGuess({0.2, 0.8});
  in->setErrMinAlpha(1.0e-4);
  in->setNIterAlpha(3);

  VSStlsManager mgr(in);
  mgr.setAlpha(0.5);
  EXPECT_DOUBLE_EQ(mgr.getAlpha(), 0.5);
  const auto &asVsManager = static_cast<const VSManager &>(mgr);
  EXPECT_FALSE(asVsManager.getWvg(GridPoints::CENTER).empty());
  EXPECT_DOUBLE_EQ(asVsManager.getCoupling(GridPoints::CENTER),
                   in->getCoupling());
  EXPECT_DOUBLE_EQ(asVsManager.getDegeneracy(GridPoints::CENTER),
                   in->getDegeneracy());
}
