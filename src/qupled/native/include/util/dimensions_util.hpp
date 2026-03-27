#ifndef DIMENSIONS_UTIL_HPP
#define DIMENSIONS_UTIL_HPP

namespace dimensionsUtil {

  enum class Dimension { D3, D2, Default = D3 };

  class DimensionsHandler {
  public:

    virtual ~DimensionsHandler() noexcept = default;
    void compute(const Dimension &dim);

  protected:

    virtual void compute2D() = 0;
    virtual void compute3D() = 0;
  };

} // namespace dimensionsUtil

#endif // DIMENSIONS_UTIL_HPP