#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <cstddef>
#include <vector>

// -----------------------------------------------------------------
// Class to represent 3D vectors (stored in row-major order)
// -----------------------------------------------------------------

class Vector3D {

public:

  // Constructors
  Vector3D(const size_t s1_, const size_t s2_, const size_t s3_)
      : v(s1_ * s2_ * s3_, 0.0),
        s1(s1_),
        s2(s2_),
        s3(s3_) {}
  explicit Vector3D()
      : Vector3D(0, 0, 0) {}

  // Return overall size of the array
  size_t size() const;

  // Return size along dimension i (0, 1 or 2)
  size_t size(const size_t i) const;

  // Check if the array is empty
  bool empty() const;

  // Resize the array to new sizes s1_, s2_, s3_
  void resize(const size_t s1_, const size_t s2_, const size_t s3_);

  // Access operator
  double &operator()(const size_t i, const size_t j, const size_t k);
  const double &
  operator()(const size_t i, const size_t j, const size_t k) const;
  const double &operator()(const size_t i, const size_t j) const;
  const double &operator()(const size_t i) const;

  // Equality operator
  bool operator==(const Vector3D &other) const;

  // Begin and end iterators
  std::vector<double>::iterator begin();
  std::vector<double>::iterator end();

  // Begin and end iterators with constant access
  std::vector<double>::const_iterator begin() const;
  std::vector<double>::const_iterator end() const;

  // Pointer to the underlying vector
  double *data();

  // Pointer to the underlying vector with constant access
  const double *data() const;

  // Set all values in the vector to num
  void fill(const double &num);
  // Set all values of row j and column k to num of a specific index i
  void fill(const size_t i, const double &num);
  // Set all entries of row i and column j to num
  void fill(const size_t i, const size_t j, const double &num);
  // Set all entries of row i and column j to be a copy of vector num
  void fill(const size_t i, const size_t j, const std::vector<double> &num);

  // Sum this vector and v_
  void sum(const Vector3D &v_);

  // Subtract this vector and v_
  void diff(const Vector3D &v_);

  // Element-wise multiplication between this vector and v_
  void mult(const Vector3D &v_);

  // Multiply each entry in the vector times num
  void mult(const double &num);

  // Element-wise division between this vector and v_
  void div(const Vector3D &v_);

  // Element-wise linear combination v += num*v_
  void linearCombination(const Vector3D &v_, const double &num_);

private:

  // Underlying vector structure with s1*s2*s3 entries
  std::vector<double> v;

  // First dimension
  size_t s1;

  // Second dimension
  size_t s2;

  // Third dimension
  size_t s3;
};

#endif
