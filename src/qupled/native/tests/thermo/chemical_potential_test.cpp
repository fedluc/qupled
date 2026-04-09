#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "schemes/input.hpp"
#include "thermo/chemical_potential.hpp"
#include "util/numerics.hpp"

TEST(ChemicalPotentialTest, StartsAsSentinelAndComputes2DClosedForm) {
  auto in = std::make_shared<Input>();
  in->setDegeneracy(1.25);

  ChemicalPotential cp(in);
  EXPECT_TRUE(std::isnan(cp.get()));

  cp.compute(dimensionsUtil::Dimension::D2);
  const double expected = std::log(std::exp(1.0 / 1.25) - 1.0);
  EXPECT_NEAR(cp.get(), expected, 1.0e-12);
}

TEST(ChemicalPotentialTest, Computes3DRootWithSmallResidual) {
  auto in = std::make_shared<Input>();
  in->setDegeneracy(2.0);
  in->setChemicalPotentialGuess({-10.0, 10.0});

  ChemicalPotential cp(in);
  cp.compute(dimensionsUtil::Dimension::D3);

  const double residual = SpecialFunctions::fermiDirac12(cp.get())
                          - 2.0 / (3.0 * std::pow(in->getDegeneracy(), 1.5));
  EXPECT_NEAR(residual, 0.0, 1.0e-8);
}
