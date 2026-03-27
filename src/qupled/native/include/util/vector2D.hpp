#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cstddef>
#include <span>
#include <vector>

/**
 * @brief Flat-storage 2D array with row-major layout.
 *
 * Stores an @f$s_1 \times s_2@f$ matrix as a contiguous @c std::vector<double>
 * in row-major (C) order.  Provides element access via @c operator()(i, j),
 * row-span access via @c operator[](i), and arithmetic helpers used throughout
 * the self-consistency iterations.
 */
class Vector2D {

public:

  /**
   * @brief Construct a zero-initialized matrix of size @p s1_ × @p s2_.
   * @param s1_ Number of rows.
   * @param s2_ Number of columns.
   */
  Vector2D(const size_t s1_, const size_t s2_)
      : v(s1_ * s2_, 0.0),
        s1(s1_),
        s2(s2_) {}

  /** @brief Construct a 0 × 0 empty array. */
  explicit Vector2D()
      : Vector2D(0, 0) {}

  /**
   * @brief Construct from a vector-of-vectors.
   * @param v_ Source 2D data; inner vectors must all have the same length.
   */
  explicit Vector2D(const std::vector<std::vector<double>> &v_);

  /**
   * @brief Construct a 1 × n array from a flat vector.
   * @param v_ Source data.
   */
  explicit Vector2D(const std::vector<double> &v_);

  /**
   * @brief Return the total number of stored elements (@f$s_1 \times s_2@f$).
   * @return Total element count.
   */
  size_t size() const;

  /**
   * @brief Return the size along dimension @p i.
   * @param i 0 for rows, 1 for columns.
   * @return Size of the requested dimension.
   */
  size_t size(const size_t i) const;

  /**
   * @brief Return true if the array has no elements.
   * @return True if @f$s_1 = 0@f$ or @f$s_2 = 0@f$.
   */
  bool empty() const;

  /**
   * @brief Resize the array, discarding existing data.
   * @param s1_ New row count.
   * @param s2_ New column count.
   */
  void resize(const size_t s1_, const size_t s2_);

  /**
   * @brief Mutable element access.
   * @param i Row index.
   * @param j Column index.
   * @return Reference to element (i, j).
   */
  double &operator()(const size_t i, const size_t j);

  /**
   * @brief Immutable element access.
   * @param i Row index.
   * @param j Column index.
   * @return Const reference to element (i, j).
   */
  const double &operator()(const size_t i, const size_t j) const;

  /**
   * @brief Return a mutable span over row @p i.
   * @param i Row index.
   * @return Span covering the @p s2 elements of row @p i.
   */
  std::span<double> operator[](const size_t i);

  /**
   * @brief Return an immutable span over row @p i.
   * @param i Row index.
   * @return Const span covering the @p s2 elements of row @p i.
   */
  std::span<const double> operator[](const size_t i) const;

  /**
   * @brief Element-wise equality comparison.
   * @param other Array to compare against.
   * @return True if all elements are equal.
   */
  bool operator==(const Vector2D &other) const;

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
   * @brief Set every element of row @p i to @p num.
   * @param i   Row index.
   * @param num Fill value.
   */
  void fill(const size_t i, const double &num);

  /**
   * @brief Set row @p i to a copy of @p num.
   * @param i   Row index.
   * @param num Source vector (must have length @p s2).
   */
  void fill(const size_t i, const std::vector<double> &num);

  /**
   * @brief Add @p v_ element-wise to this array (in-place).
   * @param v_ Array to add (must have the same shape).
   */
  void sum(const Vector2D &v_);

  /**
   * @brief Subtract @p v_ element-wise from this array (in-place).
   * @param v_ Array to subtract (must have the same shape).
   */
  void diff(const Vector2D &v_);

  /**
   * @brief Multiply element-wise by @p v_ (in-place).
   * @param v_ Array to multiply by (must have the same shape).
   */
  void mult(const Vector2D &v_);

  /**
   * @brief Multiply every element by the scalar @p num (in-place).
   * @param num Scaling factor.
   */
  void mult(const double &num);

  /**
   * @brief Divide element-wise by @p v_ (in-place).
   * @param v_ Divisor array (must have the same shape).
   */
  void div(const Vector2D &v_);

  /**
   * @brief Accumulate @f$\text{num} \cdot v\_@f$ into this array (in-place).
   * @param v_   Array to scale and add.
   * @param num_ Scaling factor.
   */
  void linearCombination(const Vector2D &v_, const double &num_);

private:

  /** @brief Flat storage for all @f$s_1 \times s_2@f$ elements (row-major). */
  std::vector<double> v;
  /** @brief Number of rows. */
  size_t s1;
  /** @brief Number of columns. */
  size_t s2;
};

#endif
