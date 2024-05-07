#ifndef BIN_UTIL_HPP
#define BIN_UTIL_HPP

#include <fstream>

// -----------------------------------------------------------------
// Utility functions to manipulate binary files
// -----------------------------------------------------------------

namespace binUtil {

  template <typename T>
  void writeNum(std::ofstream &file, const T &num) {
    file.write(reinterpret_cast<const char *>(&num), sizeof(num));
  }

  template <typename T>
  void writeDataToBinary(std::ofstream &file, const double &data) {
    writeNum<double>(file, data);
  }

  template <typename T>
  void writeDataToBinary(std::ofstream &file, const int &data) {
    writeNum<int>(file, data);
  }

  template <typename T>
  void writeDataToBinary(std::ofstream &file, const T &data) {
    for (auto &el : data) {
      writeDataToBinary<decltype(el)>(file, el);
    }
  }

  template <typename T>
  void readNum(std::ifstream &file, T &num) {
    file.read((char *)&num, sizeof(T));
  }

  template <typename T>
  void readDataFromBinary(std::ifstream &file, double &data) {
    readNum<double>(file, data);
  }

  template <typename T>
  void readDataFromBinary(std::ifstream &file, int &data) {
    readNum<int>(file, data);
  }

  template <typename T>
  void readDataFromBinary(std::ifstream &file, T &data) {
    for (auto &el : data) {
      readDataFromBinary<decltype(el)>(file, el);
    }
  }

} // namespace binUtil

#endif
