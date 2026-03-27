#ifndef FORMAT_HPP
#define FORMAT_HPP

/*
 * Disclaimer:
 * This header conditionally includes either {fmt} or std::format depending on
 * the platform. This is necessary because macOS compilers do not fully support
 * std::format introduced in C++20.
 */

#if defined(__APPLE__)
#include <fmt/core.h>
#else
#include <format>
#endif

#include <string>
#include <utility>

/**
 * @brief Platform-portable string formatting wrapper.
 *
 * On macOS, delegates to the {fmt} library because Apple's C++ standard
 * library does not ship @c std::format. On other platforms, @c std::format
 * is used directly. The calling code can use @p formatUtil::format with the
 * same syntax as @c std::format regardless of platform.
 */
namespace formatUtil {

#if defined(__APPLE__)

  /**
   * @brief Format a string using the {fmt} library.
   * @tparam Args  Argument types to format.
   * @param fmt_str Format string (must be a compile-time constant).
   * @param args    Arguments to substitute into @p fmt_str.
   * @return Formatted @c std::string.
   */
  template <typename... Args>
  std::string format(fmt::format_string<Args...> fmt_str, Args &&...args) {
    return fmt::format(fmt_str, std::forward<Args>(args)...);
  }

#else

  /**
   * @brief Format a string using @c std::format.
   * @tparam Args  Argument types to format.
   * @param fmt_str Format string (must be a compile-time constant).
   * @param args    Arguments to substitute into @p fmt_str.
   * @return Formatted @c std::string.
   */
  template <typename... Args>
  std::string format(std::format_string<Args...> fmt_str, Args &&...args) {
    return std::format(fmt_str, std::forward<Args>(args)...);
  }

#endif

} // namespace formatUtil

#endif
