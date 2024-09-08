#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cstddef>
#include <vector>

// -----------------------------------------------------------------
// Class to represent 2D vectors (stored in row-major order)
// -----------------------------------------------------------------

class Vector2D {

public:

  // Constructors
  Vector2D(const size_t s1_, const size_t s2_)
      : v(s1_ * s2_, 0.0),
        s1(s1_),
        s2(s2_) {}
  explicit Vector2D()
      : Vector2D(0, 0) {}
  explicit Vector2D(const std::vector<std::vector<double>> &v_);
  explicit Vector2D(const std::vector<double> &v_);

  // Return overall size of the array
  size_t size() const;

  // Return size along dimension i (0 or 1)
  size_t size(const size_t i) const;

  // Check if the array is empty
  bool empty() const;

  // Resize the array to new sizes s1_ and s2_
  void resize(const size_t s1_, const size_t s2_);

  // Access operator
  double &operator()(const size_t i, const size_t j);
  const double &operator()(const size_t i, const size_t j) const;
  const double &operator()(const size_t i) const;

  // Equality operator
  bool operator==(const Vector2D &other) const;

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

  // Set all entries of row i to num
  void fill(const size_t i, const double &num);

  // Set all entries of row i to be a copy of vector num
  void fill(const size_t i, const std::vector<double> &num);

  // Sum this vector and v_
  void sum(const Vector2D &v_);

  // Subtract this vector and v_
  void diff(const Vector2D &v_);

  // Element-wise multiplication between this vector and v_
  void mult(const Vector2D &v_);

  // Multiply each entry in the vector times num
  void mult(const double &num);

  // Element-wise division between this vector and v_
  void div(const Vector2D &v_);

  // Element-wise linear combination v += num*v_
  void linearCombination(const Vector2D &v_, const double &num_);

private:

  // Underlying vector structure with s1*s2 entries
  std::vector<double> v;

  // First dimension (rows)
  size_t s1;

  // Second dimension (columns)
  size_t s2;
};

#endif
