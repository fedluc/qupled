#ifndef TEST_FIXTURES_INPUT_BUILDERS_HPP
#define TEST_FIXTURES_INPUT_BUILDERS_HPP

#include <memory>
#include <string>

#include "schemes/input.hpp"

namespace testFixtures {

inline void configureBaseInput(Input &in,
                               const std::string &theory,
                               const dimensionsUtil::Dimension dim,
                               const double coupling,
                               const double degeneracy,
                               const int nMatsubara = 4) {
  in.setTheory(theory);
  in.setDimension(dim);
  in.setCoupling(coupling);
  in.setDegeneracy(degeneracy);
  in.setIntError(1.0e-6);
  in.setNMatsubara(nMatsubara);
  in.setWaveVectorGridRes(0.5);
  in.setWaveVectorGridCutoff(2.0);
  in.setFrequencyCutoff(20.0);
  in.setNThreads(1);
  in.setInt2DScheme("full");
}

inline std::shared_ptr<Input> makeBaseInput(
    const std::string &theory = "HF",
    const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
    const double coupling = 1.0,
    const double degeneracy = 0.5,
    const int nMatsubara = 4) {
  auto in = std::make_shared<Input>();
  configureBaseInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
  return in;
}

} // namespace testFixtures

#endif
