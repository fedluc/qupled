#include <numeric>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "numerics.hpp"
#include "util.hpp"

using namespace std;

namespace numUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods
  // -----------------------------------------------------------------
  
  bool equalTol(const double& x, const double &y) {
    return abs(x - y) < x * dtol;
  }

  bool largerThan(const double& x, const double &y) {
    return x - y > x * dtol;
  }
  
}

namespace vecUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods
  // -----------------------------------------------------------------
  
  // Element-wise sum between two vectors
  vector<double> sum(const vector<double> &v1,
		     const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(),
	      v2.begin(), res.begin(),
	      plus<double>());
    return res;
  }

  // Element-wise difference between two vectors
  vector<double> diff(const vector<double> &v1,
		      const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(),
	      v2.begin(), res.begin(),
	      minus<double>());
    return res;
  }

  // Element-wise multiplication between two vectors
  vector<double> mult(const vector<double> &v1,
		      const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(),
	      v2.begin(), res.begin(),
	      multiplies<double>());
    return res;
  }

  // Element-wise multiplication between two vectors
  vector<double> div(const vector<double> &v1,
		     const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(),
	      v2.begin(), res.begin(),
	      divides<double>());
    return res;
  }
  
  // Element-wise multiplication of a vector and a scalar
  vector<double> mult(const vector<double> &v,
		      const double a) {
    vector<double> res = v;
    transform(res.begin(), res.end(), res.begin(), [&a](double c){return c*a;});
    return res;
  }
  
  // Root square difference between two vectors
  double rms(const vector<double> &v1,
	     const vector<double> &v2,
	     const bool normalize) {
    const vector<double> tmp = diff(v1,v2);
    double rms = inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.0);
    if (normalize) rms /= tmp.size();
    return sqrt(rms);
  }

  // Fill vector with constant value
  void fill(vector<double> &v,
	    const double num) {
   std::for_each(v.begin(), v.end(), [&](double &vi){ vi = num;}); 
  }
  
  // -----------------------------------------------------------------
  // Vector2D class
  // -----------------------------------------------------------------
  
  Vector2D::Vector2D(const vector<vector<double>>& v_)  {
    s1 = v_.size();
    s2 = (s1 > 0) ? v_[0].size() : 0;
    v = vector<double>(s1 * s2, 0.0);
    size_t cnt = 0;
    for (const auto& vi : v_) {
      assert(vi.size() == s2);
      std::copy(vi.begin(), vi.end(), v.begin() + cnt);
      cnt += s2;
    }
  }
  
  size_t Vector2D::size() const {
    return s1*s2;
  }
  
  size_t Vector2D::size(const size_t i) const {
    assert(i == 0 || i == 1);
    return (i == 0) ? s1 : s2;
  }

  bool Vector2D::empty() const {
    return v.empty();
  }
  
  void Vector2D::resize(const size_t s1_, const size_t s2_) {
    v.clear();
    s1 = s1_;
    s2 = s2_;
    v.resize(s1_*s2_, 0.0);
  }

  double& Vector2D::operator()(const size_t i, const size_t j) {
    return v[j + i*s2];
  }
  
  const double& Vector2D::operator()(const size_t i, const size_t j) const {
    return v[j + i*s2];
  }

  const double& Vector2D::operator()(const size_t i) const {
    return operator()(i,0);
  }

  bool Vector2D::operator==(const Vector2D& other) const {
    return v == other.v && s1 == other.s1 && s2 == other.s2;
  }
  
  vector<double>::iterator Vector2D::begin() {
    return v.begin();
  }

  vector<double>::iterator Vector2D::end() {
    return v.end();
  }

  vector<double>::const_iterator Vector2D::begin() const {
    return v.begin();
  }

  vector<double>::const_iterator Vector2D::end() const {
    return v.end();
  }
  
  
  void Vector2D::fill(const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi){ vi = num;});
  }

  void Vector2D::fill(const size_t i, const double &num) {
    const auto &dest = v.begin() + i*s2;
    std::for_each(dest, dest + s2, [&](double &vi){ vi = num;});
  }
  
  void Vector2D::fill(const size_t i, const vector<double> &num) {
    assert(num.size() == s2);
    std::copy(num.begin(), num.end(), v.begin() + i*s2);
  }
  
  void Vector2D::sum(const Vector2D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::sum(v, v_.v);
  }

  void Vector2D::diff(const Vector2D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::diff(v, v_.v);
  }

  void Vector2D::mult(const Vector2D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::mult(v, v_.v);
  }

  void Vector2D::mult(const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi){ vi *= num;});
  }

  void Vector2D::div(const Vector2D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::div(v, v_.v);
  }

  // -----------------------------------------------------------------
  // Vector3D class
  // -----------------------------------------------------------------
  
  size_t Vector3D::size() const {
    return s1*s2*s3;
  } 
  
  size_t Vector3D::size(const size_t i) const {
    assert(i == 0 || i == 1 || i == 2);
    if (i == 0) return s1;
    if (i == 1) return s2;
    return s3;
  }

  bool Vector3D::empty() const {
    return v.empty();
  }

  void Vector3D::resize(const size_t s1_,
			const size_t s2_,
			const size_t s3_) {
    v.clear();
    s1 = s1_;
    s2 = s2_;
    s3 = s3_;
    v.resize(s1_*s2_*s3_, 0.0);
  }

  double& Vector3D::operator()(const size_t i,
			       const size_t j,
			       const size_t k) {
    return v[k + j*s3 + i*s2*s3];
  }
  
  const double& Vector3D::operator()(const size_t i,
				     const size_t j,
				     const size_t k) const {
    return v[k + j*s3 + i*s2*s3];
  }

  const double& Vector3D::operator()(const size_t i,
				     const size_t j) const {
    return operator()(i, j, 0);
  }

  const double& Vector3D::operator()(const size_t i) const {
    return operator()(i, 0);
  }

  bool Vector3D::operator==(const Vector3D& other) const {
    return v == other.v && s1 == other.s1
      && s2 == other.s2 && s3 == other.s3;
  }
  
  vector<double>::iterator Vector3D::begin() {
    return v.begin();
  }

  vector<double>::iterator Vector3D::end() {
    return v.end();
  }

  vector<double>::const_iterator Vector3D::begin() const {
    return v.begin();
  }

  vector<double>::const_iterator Vector3D::end() const {
    return v.end();
  }

  void Vector3D::fill(const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi){ vi = num;});
  }

  void Vector3D::fill(const size_t i,
		      const size_t j,
		      const double &num) {
    const auto &dest = v.begin() + j*s3 + i*s2*s3;
    std::for_each(dest, dest + s3, [&](double &vi){ vi = num;});
  }
  
  void Vector3D::fill(const size_t i,
		      const size_t j,
		      const vector<double> &num) {
    assert(num.size() == s3);
    std::copy(num.begin(), num.end(), v.begin() + j*s3 + i*s2*s3);
  }  

  void Vector3D::sum(const Vector3D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::sum(v, v_.v);
  }

  void Vector3D::diff(const Vector3D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::diff(v, v_.v);
  }

  void Vector3D::mult(const Vector3D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::mult(v, v_.v);
  }

  void Vector3D::mult(const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi){ vi *= num;});
  }

  void Vector3D::div(const Vector3D &v_) {
    assert(v_.size() == v.size());
    v = vecUtil::div(v, v_.v);
  }

  // -----------------------------------------------------------------
  // Stand-alone methods to convert between Python and C++ arrays
  // -----------------------------------------------------------------
  
  void python::CheckRowMajor(const bn::ndarray &nda) {
    const bn::ndarray::bitflag flags = nda.get_flags();
    const bool isRowMajor = flags & bn::ndarray::C_CONTIGUOUS;
    if (!isRowMajor) {
      MPIUtil::throwError("The numpy array is not stored in row major order (c-contiguous)");
    }
  }
  
  vector<double> python::toVector(const bn::ndarray &nda){
    if (nda.get_nd() != 1) {
      MPIUtil::throwError("Incorrect numpy array dimensions");
    }
    const Py_intptr_t* shape = nda.get_shape();
    const int dim = nda.get_nd();
    // the numpy array is flattened to a one dimensional std::vector
    Py_intptr_t n = 1;
    for (int i = 0; i < dim; ++i){ n *= shape[i]; }
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) { v[i] = *(ptr + i); }
    return v;
  }

  vector<double> python::toVector(const bp::list &list){
    int n = len(list);
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i){ v[i] = bp::extract<double>(list[i]); }
    return v;
  }

  Vector2D python::toVector2D(const bn::ndarray &nda){
    if (nda.get_nd() != 2) {
      MPIUtil::throwError("Incorrect numpy array dimensions");
    }
    CheckRowMajor(nda);
    const Py_intptr_t* shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    Vector2D v(sz1, sz2);
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    for (int i = 0; i < sz1; ++i){
      for (int j = 0; j < sz2; ++j) {
	v(i,j) = *(ptr + j + i*sz2);
      }
    }
    return v;
  }

  vector<vector<double>> python::toDoubleVector(const bn::ndarray &nda){
    if (nda.get_nd() != 2) {
      MPIUtil::throwError("Incorrect numpy array dimensions");
    }
    CheckRowMajor(nda);
    const Py_intptr_t* shape = nda.get_shape();
    const int sz1 = shape[0];
    const int sz2 = shape[1];
    vector<vector<double>> v(sz1);
    double* ptr = reinterpret_cast<double*>(nda.get_data());
    for (int i = 0; i < sz1; ++i){
      v[i].resize(sz2);
      for (int j = 0; j < sz2; ++j) {
	v[i][j] = *(ptr + j + i*sz2);
      }
    }
    return v;
  }
  
  template<typename T>
  bn::ndarray python::toNdArrayT(const T &v){
    Py_intptr_t shape[1];
    shape[0] = v.size();
    bn::ndarray result = bn::zeros(1, shape, bn::dtype::get_builtin<double>());
    std::copy(v.begin(), v.end(), reinterpret_cast<double*>(result.get_data()));
    return result;
  }

  bn::ndarray python::toNdArray(const vector<double> &v){
    return toNdArrayT(v);
  }
  
  bn::ndarray python::toNdArray2D(const Vector2D &v){
    bn::ndarray result = toNdArrayT(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1)));
    return result;
  }

  bn::ndarray python::toNdArray2D(const vector<vector<double>> &v){
    return toNdArray2D(Vector2D(v));
  }
  
  bn::ndarray python::toNdArray3D(const Vector3D &v){
    bn::ndarray result = toNdArrayT(v);
    result = result.reshape(bp::make_tuple(v.size(0), v.size(1), v.size(2)));
    return result;
  }
  
}

namespace thermoUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods
  // -----------------------------------------------------------------

  double computeInternalEnergy(const vector<double> &wvg,
			       const vector<double> &ssf,
			       const double &coupling) {
    const Interpolator1D itp(wvg, ssf);
    Integrator1D itg;
    const InternalEnergy uInt(coupling, wvg.front(), wvg.back(), itp, itg);
    return uInt.get();
  }

    
  double computeFreeEnergy(const vector<double> &grid,
			   const vector<double> &rsu,
			   const double &coupling) {
    return computeFreeEnergy(grid, rsu, coupling, true);
  }
  
  double computeFreeEnergy(const vector<double> &grid,
			   const vector<double> &rsu,
			   const double &coupling,
			   const bool normalize) {
    if (numUtil::largerThan(coupling, grid.back())) {
      MPIUtil::throwError("The coupling parameter is out of range"
			  " for the current grid, the free energy cannot be computed");
    }
    const Interpolator1D itp(grid, rsu);
    Integrator1D itg;
    const FreeEnergy freeEnergy(coupling, itp, itg, normalize);
    return freeEnergy.get();
  }
  
  vector<double> computeRdf(const vector<double> &r,
			    const vector<double> &wvg,
			    const vector<double> &ssf) {
    assert(ssf.size() > 0 && wvg.size() > 0);
    const Interpolator1D itp(wvg, ssf);
    const int nr = r.size();
    vector<double> rdf(nr);
    Integrator1DFourier itgf(0.0);
    Integrator1D itg(1.0e-6);
    for (int i=0; i<nr; ++i){
      const Rdf rdfTmp(r[i], wvg.back(), itp, itg, itgf);
      rdf[i] = rdfTmp.get();
    }
    return rdf;
  }
  
  // -----------------------------------------------------------------
  // InternalEnergy class
  // -----------------------------------------------------------------
  
  double InternalEnergy::ssf(double y) const {
    return ssfi.eval(y);
  }
  
  double InternalEnergy::integrand(double y) const {
    return ssf(y) - 1;
  }

  double InternalEnergy::get() const  {
    auto func = [&](double y)->double{return integrand(y);};
    itg.compute(func, yMin, yMax);
    return itg.getSolution()/(M_PI * rs * lambda);
  }

  // -----------------------------------------------------------------
  // FreeEnergy class
  // -----------------------------------------------------------------
  
  double FreeEnergy::get() const  {
    auto func = [&](double y)->double{return rsui.eval(y);};
    itg.compute(func, 0.0, rs);
    if (normalize) { return (rs == 0.0) ? -numUtil::Inf : itg.getSolution()/rs/rs; };
    return itg.getSolution();
  }
  
  // -----------------------------------------------------------------
  // Rdf class
  // -----------------------------------------------------------------
  
  double Rdf::ssf(double y) const {
    return ssfi.eval(y);
  }

  double Rdf::integrand(double y) const {
    if (y > cutoff) return 0;
    const double yssf = y * (ssf(y) - 1);
    return (r == 0.0) ? y * yssf : yssf;
  }

  double Rdf::get() const {
    auto func = [&](double y)->double{return integrand(y);};
    if (r == 0) {
      itg.compute(func, 0.0, cutoff);
      return 1 + 1.5 * itg.getSolution();
    }
    else {
      itgf.setR(r);
      itgf.compute(func);
      return 1 + 1.5 * itgf.getSolution()/r;
    }
  }
  
}
  
namespace MPIUtil {

  // -----------------------------------------------------------------
  // Stand-alone methods
  // -----------------------------------------------------------------

  int getRank() {
    int rank;
    MPI_Comm_rank(MPICommunicator, &rank);
    return rank;
  }

  bool isRoot() {
    return getRank() == 0;
  }

  bool isSingleProcess() {
    int numRanks;
    MPI_Comm_size(MPICommunicator, &numRanks);
    return numRanks == 1;
  }

  void throwError(const string& errMsg) {
    if (MPIUtil::isSingleProcess()) {
      // Throw a catchable error if only one process is used
      throw std::runtime_error(errMsg);
    }
    // Abort MPI if more than one process is running
    std::cerr << errMsg << std::endl;
    MPIUtil::abort();
  }
  
  void abort() {
    MPI_Abort(MPICommunicator, 1);
  }

  double timer() {
    return MPI_Wtime();
  }
  
}
