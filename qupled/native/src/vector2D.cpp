#include "vector2D.hpp"
#include "vector_util.hpp"
#include <algorithm>
#include <cassert>

using namespace std;

Vector2D::Vector2D(const vector<vector<double>> &v_) {
  s1 = v_.size();
  s2 = (s1 > 0) ? v_[0].size() : 0;
  v = vector<double>(s1 * s2, 0.0);
  size_t cnt = 0;
  for (const auto &vi : v_) {
    assert(vi.size() == s2);
    std::copy(vi.begin(), vi.end(), v.begin() + cnt);
    cnt += s2;
  }
}

Vector2D::Vector2D(const vector<double> &v_) {
  s1 = v_.size();
  s2 = 1;
  v = v_;
}

size_t Vector2D::size() const { return s1 * s2; }

size_t Vector2D::size(const size_t i) const {
  assert(i == 0 || i == 1);
  return (i == 0) ? s1 : s2;
}

bool Vector2D::empty() const { return v.empty(); }

void Vector2D::resize(const size_t s1_, const size_t s2_) {
  v.clear();
  s1 = s1_;
  s2 = s2_;
  v.resize(s1_ * s2_, 0.0);
}

double &Vector2D::operator()(const size_t i, const size_t j) {
  return v[j + i * s2];
}

const double &Vector2D::operator()(const size_t i, const size_t j) const {
  return v[j + i * s2];
}

const double &Vector2D::operator()(const size_t i) const {
  return operator()(i, 0);
}

bool Vector2D::operator==(const Vector2D &other) const {
  return v == other.v && s1 == other.s1 && s2 == other.s2;
}

vector<double>::iterator Vector2D::begin() { return v.begin(); }

vector<double>::iterator Vector2D::end() { return v.end(); }

vector<double>::const_iterator Vector2D::begin() const { return v.begin(); }

vector<double>::const_iterator Vector2D::end() const { return v.end(); }

double *Vector2D::data() { return v.data(); }

const double *Vector2D::data() const { return v.data(); }

void Vector2D::fill(const double &num) {
  std::for_each(v.begin(), v.end(), [&](double &vi) { vi = num; });
}

void Vector2D::fill(const size_t i, const double &num) {
  const auto &dest = v.begin() + i * s2;
  std::for_each(dest, dest + s2, [&](double &vi) { vi = num; });
}

void Vector2D::fill(const size_t i, const vector<double> &num) {
  assert(num.size() == s2);
  std::copy(num.begin(), num.end(), v.begin() + i * s2);
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
  std::for_each(v.begin(), v.end(), [&](double &vi) { vi *= num; });
}

void Vector2D::div(const Vector2D &v_) {
  assert(v_.size() == v.size());
  v = vecUtil::div(v, v_.v);
}

void Vector2D::linearCombination(const Vector2D &v_, const double &num) {
  assert(v_.size() == v.size());
  v = vecUtil::linearCombination(v, 1, v_.v, num);
}
