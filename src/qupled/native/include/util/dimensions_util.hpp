#ifndef DIMENSIONS_UTIL_HPP
#define DIMENSIONS_UTIL_HPP

/**
 * @brief Dimension-aware dispatch infrastructure for 2D and 3D calculations.
 *
 * Provides an enumeration of supported spatial dimensions and a base class
 * whose @p compute() method dispatches to dimension-specific virtual overrides.
 */
namespace dimensionsUtil {

  /** @brief Supported spatial dimensions. @p Default maps to 3D. */
  enum class Dimension { D3, D2, Default = D3 };

  /**
   * @brief Base class that dispatches computations based on the spatial
   * dimension.
   *
   * Derived classes implement @p compute2D() and @p compute3D() to provide
   * dimension-specific logic.  Call @p compute(dim) to invoke the correct
   * override at runtime.
   */
  class DimensionsHandler {
  public:

    /** @brief Virtual destructor. */
    virtual ~DimensionsHandler() noexcept = default;

    /**
     * @brief Dispatch to the appropriate dimension-specific compute method.
     * @param dim Spatial dimension to use for the computation.
     */
    void compute(const Dimension &dim);

  protected:

    /** @brief Compute for a 2D system. Implemented by derived classes. */
    virtual void compute2D() = 0;

    /** @brief Compute for a 3D system. Implemented by derived classes. */
    virtual void compute3D() = 0;
  };

} // namespace dimensionsUtil

#endif // DIMENSIONS_UTIL_HPP
