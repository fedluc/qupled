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
    in.setChemicalPotentialGuess({-10.0, 10.0});
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

  inline void configureIterationInput(IterationInput &in,
                                      const std::string &theory,
                                      const dimensionsUtil::Dimension dim,
                                      const double coupling,
                                      const double degeneracy,
                                      const int nMatsubara = 2) {
    configureBaseInput(in, theory, dim, coupling, degeneracy, nMatsubara);
    in.setErrMin(1.0e-4);
    in.setNIter(2);
    in.setMixingParameter(0.5);
  }

  inline std::shared_ptr<StlsInput> makeStlsInput(
      const std::string &theory = "STLS",
      const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
      const double coupling = 1.0,
      const double degeneracy = 0.7,
      const int nMatsubara = 2) {
    auto in = std::make_shared<StlsInput>();
    configureIterationInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
    return in;
  }

  inline std::shared_ptr<StlsIetInput> makeStlsIetInput(
      const std::string &theory = "STLS-HNC",
      const std::string &mapping = "sqrt",
      const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
      const double coupling = 1.0,
      const double degeneracy = 0.7,
      const int nMatsubara = 2) {
    auto in = std::make_shared<StlsIetInput>();
    configureIterationInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
    in->setMapping(mapping);
    return in;
  }

  inline std::shared_ptr<QstlsInput> makeQstlsInput(
      const std::string &theory = "QSTLS",
      const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
      const double coupling = 1.0,
      const double degeneracy = 0.7,
      const int nMatsubara = 2) {
    auto in = std::make_shared<QstlsInput>();
    configureIterationInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
    return in;
  }

  inline std::shared_ptr<QstlsIetInput> makeQstlsIetInput(
      const std::string &theory = "QSTLS-HNC",
      const std::string &mapping = "sqrt",
      const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
      const double coupling = 1.0,
      const double degeneracy = 0.7,
      const int nMatsubara = 2) {
    auto in = std::make_shared<QstlsIetInput>();
    configureIterationInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
    in->setMapping(mapping);
    return in;
  }

  inline std::shared_ptr<QVSStlsInput> makeQVSStlsInput(
      const std::string &theory = "QVSSTLS",
      const dimensionsUtil::Dimension dim = dimensionsUtil::Dimension::D3,
      const double coupling = 1.0,
      const double degeneracy = 0.7,
      const int nMatsubara = 2) {
    auto in = std::make_shared<QVSStlsInput>();
    configureIterationInput(*in, theory, dim, coupling, degeneracy, nMatsubara);
    in->setCouplingResolution(0.5);
    in->setDegeneracyResolution(0.2);
    in->setAlphaGuess({0.2, 0.8});
    in->setErrMinAlpha(1.0e-4);
    in->setNIterAlpha(5);
    return in;
  }

} // namespace testFixtures

#endif
