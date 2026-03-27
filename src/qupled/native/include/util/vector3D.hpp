#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <cstddef>
#include <span>
#include <vector>

/**
 * @brief Flat-storage 3D array with row-major layout.
 *
 * Stores an @f$s_1 \times s_2 \times s_3@f$ array as a contiguous
 * @c std::vector<double> in row-major (C) order.  Elements are accessed via
 * @c operator()(i, j, k).  Used primarily to store the fixed component of the
 * qSTLS auxiliary density response.
 */
class Vector3D {

public:

  /**
   * @brief Construct a zero-initialized array of size @p s1_ × @p s2_ × @p s3_.
   * @param s1_ Size of the first dimension.
   * @param s2_ Size of the second dimension.
   * @param s3_ Size of the third dimension.
   */
  Vector3D(const size_t s1_, const size_t s2_, const size_t s3_)
      : v(s1_ * s2_ * s3_, 0.0),
        s1(s1_),
        s2(s2_),
        s3(s3_) {}

  /** @brief Construct a 0 × 0 × 0 empty array. */
  explicit Vector3D()
      : Vector3D(0, 0, 0) {}

  /**
   * @brief Return the total number of stored elements.
   * @return @f$s_1 \times s_2 \times s_3@f$.
   */
  size_t size() const;

  /**
   * @brief Return the size along dimension @p i.
   * @param i 0, 1, or 2 for the first, second, or third dimension.
   * @return Size of the requested dimension.
   */
  size_t size(const size_t i) const;

  /**
   * @brief Return true if the array has no elements.
   * @return True if any dimension size is zero.
   */
  bool empty() const;

  /**
   * @brief Resize the array, discarding existing data.
   * @param s1_ New size of the first dimension.
   * @param s2_ New size of the second dimension.
   * @param s3_ New size of the third dimension.
   */
  void resize(const size_t s1_, const size_t s2_, const size_t s3_);

  /**
   * @brief Mutable element access.
   * @param i First-dimension index.
   * @param j Second-dimension index.
   * @param k Third-dimension index.
   * @return Reference to element (i, j, k).
   */
  double &operator()(const size_t i, const size_t j, const size_t k);

  /**
   * @brief Immutable element access.
   * @param i First-dimension index.
   * @param j Second-dimension index.
   * @param k Third-dimension index.
   * @return Const reference to element (i, j, k).
   */
  const double &
  operator()(const size_t i, const size_t j, const size_t k) const;

  /**
   * @brief Element-wise equality comparison.
   * @param other Array to compare against.
   * @return True if all elements are equal.
   */
  bool operator==(const Vector3D &other) const;

  /** @brief Return an iterator to the first element. */
  std::vector<double>::iterator begin();
  /** @brief Return an iterator past the last element. */
  std::vector<double>::iterator end();
  /** @brief Return a const iterator to the first element. */
  std::vector<double>::const_iterator begin() const;
  /** @brief Return a const iterator past the last element. */
  std::vector<double>::const_iterator end() const;

  /**
   * @brief Return a mutable pointer to the underlying contiguous storage.
   * @return Pointer to the first element.
   */
  double *data();

  /**
   * @brief Return an immutable pointer to the underlying contiguous storage.
   * @return Const pointer to the first element.
   */
  const double *data() const;

  /**
   * @brief Set every element to @p num.
   * @param num Fill value.
   */
  void fill(const double &num);

  /**
   * @brief Set every element with first index @p i to @p num.
   * @param i   First-dimension index.
   * @param num Fill value.
   */
  void fill(const size_t i, const double &num);

  /**
   * @brief Set every element with indices (@p i, @p j, *) to @p num.
   * @param i   First-dimension index.
   * @param j   Second-dimension index.
   * @param num Fill value.
   */
  void fill(const size_t i, const size_t j, const double &num);

  /**
   * @brief Set the row (i, j, *) to a copy of @p num.
   * @param i   First-dimension index.
   * @param j   Second-dimension index.
   * @param num Source vector (must have length @p s3).
   */
  void fill(const size_t i, const size_t j, const std::vector<double> &num);

  /**
   * @brief Add @p v_ element-wise to this array (in-place).
   * @param v_ Array to add (must have the same shape).
   */
  void sum(const Vector3D &v_);

  /**
   * @brief Subtract @p v_ element-wise from this array (in-place).
   * @param v_ Array to subtract (must have the same shape).
   */
  void diff(const Vector3D &v_);

  /**
   * @brief Multiply element-wise by @p v_ (in-place).
   * @param v_ Array to multiply by (must have the same shape).
   */
  void mult(const Vector3D &v_);

  /**
   * @brief Multiply every element by the scalar @p num (in-place).
   * @param num Scaling factor.
   */
  void mult(const double &num);

  /**
   * @brief Divide element-wise by @p v_ (in-place).
   * @param v_ Divisor array (must have the same shape).
   */
  void div(const Vector3D &v_);

  /**
   * @brief Accumulate @f$\text{num} \cdot v\_@f$ into this array (in-place).
   * @param v_   Array to scale and add.
   * @param num_ Scaling factor.
   */
  void linearCombination(const Vector3D &v_, const double &num_);

private:

  /** @brief Flat storage for all @f$s_1 \times s_2 \times s_3@f$ elements
   * (row-major). */
  std::vector<double> v;
  /** @brief Size of the first dimension. */
  size_t s1;
  /** @brief Size of the second dimension. */
  size_t s2;
  /** @brief Size of the third dimension. */
  size_t s3;
};

#endif
