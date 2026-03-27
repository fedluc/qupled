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

namespace formatUtil {

#if defined(__APPLE__)

  template <typename... Args>
  std::string format(fmt::format_string<Args...> fmt_str, Args &&...args) {
    return fmt::format(fmt_str, std::forward<Args>(args)...);
  }

#else

  template <typename... Args>
  std::string format(std::format_string<Args...> fmt_str, Args &&...args) {
    return std::format(fmt_str, std::forward<Args>(args)...);
  }

#endif

} // namespace formatUtil

#endif