#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "python_wrappers.hpp"
#include "input.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
namespace vp = vecUtil::python;

// Methods that need wrapping to pass arrays between native and python

namespace VSStlsInputWrapper {

  bn::ndarray getAlphaGuess(VSStlsInput &in){
    return vp::toNdArray(in.getAlphaGuess());
  }
  
  void setAlphaGuess(VSStlsInput &in,
		     const bp::list &alphaGuess){
    in.setAlphaGuess(vp::toVector(alphaGuess));
  }
  
  struct FreeEnergyIntegrand {
    bn::ndarray grid = vp::toNdArray(vector<double>(0));
    bn::ndarray integrand = vp::toNdArray(vector<double>(0));
  };
  
  VSStlsInputWrapper::FreeEnergyIntegrand getFreeEnergyIntegrand(VSStlsInput &in){
    VSStlsInput::FreeEnergyIntegrand fxcIntegrand_ = in.getFreeEnergyIntegrand();
    VSStlsInputWrapper::FreeEnergyIntegrand fxcIntegrand;
    fxcIntegrand.grid = vp::toNdArray(fxcIntegrand_.grid);
    fxcIntegrand.integrand = vp::toNdArray2D(fxcIntegrand_.integrand);
    return fxcIntegrand;
  }
  
  void setFreeEnergyIntegrand(VSStlsInput &in,
			      const VSStlsInputWrapper::FreeEnergyIntegrand &fxcIntegrand){
    VSStlsInput::FreeEnergyIntegrand fxcIntegrand_;
    fxcIntegrand_.grid = vp::toVector(fxcIntegrand.grid);
    fxcIntegrand_.integrand = vp::toDoubleVector(fxcIntegrand.integrand);
    in.setFreeEnergyIntegrand(fxcIntegrand_);
  }
  
}

namespace QstlsInputWrapper {
  
  struct QstlsGuess {
    bn::ndarray wvg = vp::toNdArray(vector<double>(0));
    bn::ndarray ssf = vp::toNdArray(vector<double>(0));
    bn::ndarray adr = vp::toNdArray(vector<double>(0));
    int matsubara = 0;
  };

  QstlsInputWrapper::QstlsGuess getGuess(QstlsInput &in){
    QstlsInput::QstlsGuess guess_ = in.getGuess();
    QstlsInputWrapper::QstlsGuess guess;
    guess.wvg = vp::toNdArray(guess_.wvg);
    guess.ssf = vp::toNdArray(guess_.ssf);
    bn::ndarray adrTmp = vp::toNdArray2D(guess_.adr);
    guess.matsubara = guess_.matsubara;
    return guess;
  }
  
  void setGuess(QstlsInput &in,
		const QstlsInputWrapper::QstlsGuess &guess){
    QstlsInput::QstlsGuess guess_;
    guess_.wvg = vp::toVector(guess.wvg);
    guess_.ssf = vp::toVector(guess.ssf);
    if (guess.adr.shape(0) > 0) {
      guess_.adr = vp::toVector2D(guess.adr);
    }
    guess_.matsubara = guess.matsubara;
    in.setGuess(guess_);
  }
  
}

// Classes exposed to Python
BOOST_PYTHON_MODULE(qupled)
{

  // Docstring formatting
  bp::docstring_options docopt;
  docopt.enable_all();
  docopt.disable_cpp_signatures();
	
  // Numpy library initialization
  bn::initialize();

  // Class for the input of the Rpa scheme
  bp::class_<RpaInput>("RpaInput")
    .add_property("coupling",
		  &RpaInput::getCoupling,
		  &RpaInput::setCoupling)
    .add_property("degeneracy",
		  &RpaInput::getDegeneracy,
		  &RpaInput::setDegeneracy)
    .add_property("int2DScheme",
		  &RpaInput::getInt2DScheme,
		  &RpaInput::setInt2DScheme)
    .add_property("intError",
		  &RpaInput::getIntError,
		  &RpaInput::setIntError)
    .add_property("threads",
		  &RpaInput::getNThreads,
		  &RpaInput::setNThreads)
    .add_property("theory",
		  &RpaInput::getTheory,
		  &RpaInput::setTheory)
    .add_property("chemicalPotential",
		  &PyInput::getChemicalPotentialGuess,
		  &PyInput::setChemicalPotentialGuess)
    .add_property("matsubara",
		  &RpaInput::getNMatsubara,
		  &RpaInput::setNMatsubara)
    .add_property("resolution",
		  &RpaInput::getWaveVectorGridRes,
		  &RpaInput::setWaveVectorGridRes)
    .add_property("cutoff",
		  &RpaInput::getWaveVectorGridCutoff,
		  &RpaInput::setWaveVectorGridCutoff)
    .def("print", &RpaInput::print)
    .def("isEqual", &RpaInput::isEqual);

  // Class for the initial guess of the Stls scheme
  bp::class_<StlsInput::SlfcGuess>("SlfcGuess")
    .add_property("wvg",
		  &PyInput::getWvg,
		  &PyInput::setWvg)
    .add_property("slfc",
		  &PyInput::getSlfc,
		  &PyInput::setSlfc);

  // Class for the input of the Stls scheme
  bp::class_<StlsInput, bp::bases<RpaInput>>("StlsInput")
    .add_property("error",
		  &StlsInput::getErrMin,
		  &StlsInput::setErrMin)
    .add_property("mixing",
		  &StlsInput::getMixingParameter,
		  &StlsInput::setMixingParameter)
    .add_property("iet",
		  &StlsInput::getIETMapping,
		  &StlsInput::setIETMapping)
    .add_property("iterations",
		  &StlsInput::getNIter,
		  &StlsInput::setNIter)
    .add_property("outputFrequency",
		  &StlsInput::getOutIter,
		  &StlsInput::setOutIter)
    .add_property("recoveryFile",
		  &StlsInput::getRecoveryFileName,
		  &StlsInput::setRecoveryFileName)
    .add_property("guess",
		  &StlsInput::getGuess,
		  &StlsInput::setGuess)
    .def("print", &StlsInput::print)
    .def("isEqual", &StlsInput::isEqual);

  // Class for the free energy integrand of the VSStls scheme
  bp::class_<VSStlsInputWrapper::FreeEnergyIntegrand>("FreeEnergyIntegrand")
    .def_readwrite("grid", &VSStlsInputWrapper::FreeEnergyIntegrand::grid)
    .def_readwrite("integrand", &VSStlsInputWrapper::FreeEnergyIntegrand::integrand);

  // Class for the input of the VSStls scheme
  bp::class_<VSStlsInput, bp::bases<StlsInput>>("VSStlsInput")
    .add_property("errorAlpha",
		  &VSStlsInput::getErrMinAlpha,
		  &VSStlsInput::setErrMinAlpha)
    .add_property("iterationsAlpha",
		  &VSStlsInput::getNIterAlpha,
		  &VSStlsInput::setNIterAlpha)
    .add_property("alpha",
		  VSStlsInputWrapper::getAlphaGuess,
		  VSStlsInputWrapper::setAlphaGuess)
    .add_property("couplingResolution",
		  &VSStlsInput::getCouplingResolution,
		  &VSStlsInput::setCouplingResolution)
    .add_property("degeneracyResolution",
		  &VSStlsInput::getDegeneracyResolution,
		  &VSStlsInput::setDegeneracyResolution)
    .add_property("freeEnergyIntegrand",
		  VSStlsInputWrapper::getFreeEnergyIntegrand,
		  VSStlsInputWrapper::setFreeEnergyIntegrand)
    .def("print", &VSStlsInput::print)
    .def("isEqual", &VSStlsInput::isEqual);

  // Class for the initial guess of the Qstls scheme
  bp::class_<QstlsInputWrapper::QstlsGuess>("QstlsGuess")
    .def_readwrite("wvg", &QstlsInputWrapper::QstlsGuess::wvg)
    .def_readwrite("ssf", &QstlsInputWrapper::QstlsGuess::ssf)
    .def_readwrite("adr", &QstlsInputWrapper::QstlsGuess::adr)
    .def_readwrite("matsubara", &QstlsInputWrapper::QstlsGuess::matsubara);

  // Class for the input of the Qstls scheme
  bp::class_<QstlsInput, bp::bases<StlsInput>>("QstlsInput")
    .add_property("guess",
		  QstlsInputWrapper::getGuess,
		  QstlsInputWrapper::setGuess)
    .add_property("fixed",
		  &QstlsInput::getFixed,
		  &QstlsInput::setFixed)
    .add_property("fixediet",
		  &QstlsInput::getFixedIet,
		  &QstlsInput::setFixedIet)
    .def("print", &QstlsInput::print)
    .def("isEqual", &QstlsInput::isEqual);

  // Class to solve the classical RPA scheme
  bp::class_<Rpa>("Rpa",
		  bp::init<const RpaInput>())
    .def("rdf", &PyRpa::getRdf)
    .add_property("idr", &PyRpa::getIdr)
    .add_property("sdr", &PyRpa::getSdr)
    .add_property("slfc", &PyRpa::getSlfc)
    .add_property("ssf", &PyRpa::getSsf)
    .add_property("ssfHF", &PyRpa::getSsfHF)
    .add_property("uInt", &PyRpa::getUInt)
    .add_property("wvg", &PyRpa::getWvg)
    .add_property("recovery", &PyRpa::getRecoveryFileName);

  // Class to solve the classical ESA scheme
  bp::class_<ESA, bp::bases<Rpa>>("ESA",
				  bp::init<const RpaInput>());

  // Class to solve classical schemes
  bp::class_<Stls, bp::bases<Rpa>>("Stls",
				   bp::init<const StlsInput>())
    .def("compute", &PyStls::compute)
    .add_property("bf", &PyStls::getBf);

  // Class to solve the classical vs scheme
  bp::class_<VSStls, bp::bases<Rpa>>("VSStls",
				     bp::init<const VSStlsInput>())
    .def("compute", &PyVSStls::compute)
    .add_property("freeEnergyIntegrand", &PyVSStls::getFreeEnergyIntegrand)
    .add_property("freeEnergyGrid", &PyVSStls::getFreeEnergyGrid);
      
  // Class to solve quantum schemes
  bp::class_<Qstls, bp::bases<Stls>>("Qstls",
				     bp::init<const QstlsInput>())
    .def("compute", &PyQstls::compute)
    .add_property("adr", &PyQstls::getAdr);

  // Post-process methods
  bp::def("computeRdf", &PyThermo::computeRdf);
  bp::def("computeInternalEnergy", &PyThermo::computeInternalEnergy);
  bp::def("computeFreeEnergy", &PyThermo::computeFreeEnergy);
  
}
