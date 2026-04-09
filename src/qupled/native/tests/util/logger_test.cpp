#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include "util/logger.hpp"

namespace {
class TestLogger : public Logger {
public:
  explicit TestLogger(bool verbose)
      : Logger(verbose) {}

  void log(const std::string &msg) const { print(msg); }
  void logLine(const std::string &msg) const { println(msg); }
};
} // namespace

TEST(LoggerTest, PrintAndPrintlnWriteToStdoutWhenVerbose) {
  TestLogger logger(true);
  std::stringstream out;
  auto *old = std::cout.rdbuf(out.rdbuf());

  logger.log("abc");
  logger.logLine("def");

  std::cout.rdbuf(old);
  EXPECT_EQ(out.str(), "abcdef\n");
}

TEST(LoggerTest, LoggingIsMutedWhenNotVerbose) {
  TestLogger logger(false);
  std::stringstream out;
  auto *old = std::cout.rdbuf(out.rdbuf());

  logger.log("abc");
  logger.logLine("def");

  std::cout.rdbuf(old);
  EXPECT_TRUE(out.str().empty());
}
