#ifndef DIMENSIONS_UTIL_HPP
#define DIMENSIONS_UTIL_HPP

#include "input.hpp"
#include <functional>
#include <memory>

namespace dimensionsUtil {

template <typename Derived>
class DimensionAware {
protected:
    Derived& derived() { return *static_cast<Derived*>(this); }
    const Derived& derived() const { return *static_cast<const Derived*>(this); }
    
public:
void compute(Input::Dimension dim) {
    switch (dim) {
        case Input::Dimension::D2: 
            derived().compute2D(); 
            break;
        case Input::Dimension::D3: 
            derived().compute3D(); 
            break;
        default: 
            throw std::invalid_argument("Unsupported dimension");
    }
}

    virtual ~DimensionAware() noexcept = default;   

private:
    void check_dimensions(Input::Dimension dim,
                        const std::function<void()>& func2D,
                        const std::function<void()>& func3D) {
        switch (dim) {
            case Input::Dimension::D2: func2D(); break;
            case Input::Dimension::D3: func3D(); break;
            default: throw std::invalid_argument("Unsupported dimension");
        }
    }
};

} // namespace dimensionsUtil

#endif // DIMENSIONS_UTIL_HPP