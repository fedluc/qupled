#include "vector3D.hpp"
#include "vector_util.hpp"
#include <algorithm>
#include <cassert>

using namespace std;

size_t Vector3D::size() const { return s1 * s2 * s3; }

size_t Vector3D::size(const size_t i) const {
  assert(i == 0 || i == 1 || i == 2);
  if (i == 0) return s1;
  if (i == 1) return s2;
  return s3;
}

bool Vector3D::empty() const { return v.empty(); }

void Vector3D::resize(const size_t s1_, const size_t s2_, const size_t s3_) {
  v.clear();
  s1 = s1_;
  s2 = s2_;
  s3 = s3_;
  v.resize(s1_ * s2_ * s3_, 0.0);
}

double &Vector3D::operator()(const size_t i, const size_t j, const size_t k) {
  return v[k + j * s3 + i * s2 * s3];
}

const double &
Vector3D::operator()(const size_t i, const size_t j, const size_t k) const {
  return v[k + j * s3 + i * s2 * s3];
}

const double &Vector3D::operator()(const size_t i, const size_t j) const {
  return operator()(i, j, 0);
}

const double &Vector3D::operator()(const size_t i) const {
  return operator()(i, 0);
}

bool Vector3D::operator==(const Vector3D &other) const {
  return v == other.v && s1 == other.s1 && s2 == other.s2 && s3 == other.s3;
}

vector<double>::iterator Vector3D::begin() { return v.begin(); }

vector<double>::iterator Vector3D::end() { return v.end(); }

vector<double>::const_iterator Vector3D::begin() const { return v.begin(); }

vector<double>::const_iterator Vector3D::end() const { return v.end(); }

double *Vector3D::data() { return v.data(); }

const double *Vector3D::data() const { return v.data(); }

void Vector3D::fill(const double &num) {
  std::for_each(v.begin(), v.end(), [&](double &vi) { vi = num; });
}

void Vector3D::fill(const size_t i, const double &num) {
  const auto &dest = v.begin() + i * s2 * s3;
  std::for_each(dest, dest + s2 * s3, [&](double &vi) { vi = num; });
}

void Vector3D::fill(const size_t i, const size_t j, const double &num) {
  const auto &dest = v.begin() + j * s3 + i * s2 * s3;
  std::for_each(dest, dest + s3, [&](double &vi) { vi = num; });
}

void Vector3D::fill(const size_t i, const size_t j, const vector<double> &num) {
  assert(num.size() == s3);
  std::copy(num.begin(), num.end(), v.begin() + j * s3 + i * s2 * s3);
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
  std::for_each(v.begin(), v.end(), [&](double &vi) { vi *= num; });
}

void Vector3D::div(const Vector3D &v_) {
  assert(v_.size() == v.size());
  v = vecUtil::div(v, v_.v);
}

void Vector3D::linearCombination(const Vector3D &v_, const double &num) {
  assert(v_.size() == v.size());
  v = vecUtil::linearCombination(v, 1.0, v_.v, num);
}
