#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <vector>

// -----------------------------------------------------------------
// Class to represent 3D vectors (stored in row-major order)
// -----------------------------------------------------------------

class Vector3D {

public:
  
  Vector3D(const size_t s1_, const size_t s2_, const size_t s3_)
    : v(s1_ * s2_ * s3_, 0.0),
      s1(s1_),
      s2(s2_),
      s3(s3_) {}
  explicit Vector3D()
    : Vector3D(0, 0, 0) {}
  size_t size() const;
  size_t size(const size_t i) const;
  bool empty() const;
  void resize(const size_t s1_, const size_t s2_, const size_t s3_);
  double &operator()(const size_t i, const size_t j, const size_t k);
  const double &
  operator()(const size_t i, const size_t j, const size_t k) const;
  const double &operator()(const size_t i, const size_t j) const;
  const double &operator()(const size_t i) const;
  bool operator==(const Vector3D &other) const;
  std::vector<double>::iterator begin();
  std::vector<double>::iterator end();
  std::vector<double>::const_iterator begin() const;
  std::vector<double>::const_iterator end() const;
  double *data();
  const double *data() const;
  void fill(const double &num);
  void fill(const size_t i, const size_t j, const double &num);
  void fill(const size_t i, const size_t j, const std::vector<double> &num);
  void sum(const Vector3D &v_);
  void diff(const Vector3D &v_);
  void mult(const Vector3D &v_);
  void mult(const double &num);
  void div(const Vector3D &v_);
  
private:
  
  std::vector<double> v;
  size_t s1;
  size_t s2;
  size_t s3;
};

#endif
