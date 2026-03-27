#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <string>

/**
 * @brief Provides controlled stdout output with optional verbosity.
 *
 * Base class for components that need to print progress or diagnostic
 * messages. Output can be suppressed by constructing with @p verbose_ = false.
 */
class Logger {

protected:

  /**
   * @brief Construct with explicit verbosity flag.
   * @param verbose_ If false, all print calls are silenced.
   */
  Logger(const bool &verbose_)
      : verbose(verbose_) {}

  /** @brief Construct with verbosity enabled. */
  Logger()
      : Logger(true) {}

  /**
   * @brief Print a message without a trailing newline.
   * @param msg Message to print.
   */
  void print(const std::string &msg) const;

  /**
   * @brief Print a message followed by a newline.
   * @param msg Message to print.
   */
  void println(const std::string &msg) const;

private:

  /** @brief Whether messages should be printed. */
  const bool verbose;

  /** @brief Returns true if the message should be logged. */
  bool log_message() const;
};

#endif
