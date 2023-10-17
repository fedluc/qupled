#include <numeric>
#include "util.hpp"

namespace vecUtil {

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

  // Vector2D class
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

  // Vector3D class
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
  
}

namespace thermoUtil {

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

  double Rdf::ssf(double y) const {
    return ssfi.eval(y);
  }

  double Rdf::integrand(double y) const {
    if (y > cutoff) return 0;
    return y*(ssf(y) - 1);
  }

  double Rdf::get() const {
    if (r == 0) { return 0.0; }
    auto func = [&](double y)->double{return integrand(y);};
    itg.setR(r);
    itg.compute(func);
    return 1 + 1.5 * itg.getSolution()/r;
  }

  vector<double> computeRdf(const vector<double> &r,
			    const vector<double> &wvg,
			    const vector<double> &ssf) {
    assert(ssf.size() > 0 && wvg.size() > 0);
    const Interpolator1D itp(wvg, ssf);
    const int nr = r.size();
    vector<double> rdf(nr);
    Integrator1DFourier itg(0.0);
    for (int i=0; i<nr; ++i){
      const Rdf rdfTmp(r[i], wvg.back(), itp, itg);
      rdf[i] = rdfTmp.get();
    }
    return rdf;
  }

  double computeInternalEnergy(const vector<double> &wvg,
			       const vector<double> &ssf,
			       const double &coupling) {
    const Interpolator1D itp(wvg, ssf);
    Integrator1D itg;
    const InternalEnergy uInt(coupling, wvg.front(), wvg.back(), itp, itg);
    return uInt.get();
  }
  
}
  

