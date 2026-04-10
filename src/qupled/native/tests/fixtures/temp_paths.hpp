#ifndef TEST_FIXTURES_TEMP_PATHS_HPP
#define TEST_FIXTURES_TEMP_PATHS_HPP

#include <chrono>
#include <filesystem>
#include <string>

namespace testFixtures {

  inline std::filesystem::path makeUniqueTempDir(const std::string &prefix) {
    const auto stamp =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto path = std::filesystem::temp_directory_path()
                / (prefix + std::to_string(stamp));
    std::filesystem::create_directories(path);
    return path;
  }

  class ScopedTempDir {
  public:

    explicit ScopedTempDir(const std::string &prefix)
        : dir(makeUniqueTempDir(prefix)) {}

    ~ScopedTempDir() { std::filesystem::remove_all(dir); }

    const std::filesystem::path &path() const { return dir; }

  private:

    std::filesystem::path dir;
  };

} // namespace testFixtures

#endif
