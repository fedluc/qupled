#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <string>

// -----------------------------------------------------------------
// Logger class used to print information on screen
// -----------------------------------------------------------------

class Logger {

protected:

  // Constructor
  Logger(const bool &verbose_)
      : verbose(verbose_) {}
  // Output verbosity
  const bool verbose;
  // Print info on screen
  void print(const std::string &msg) const;
  void println(const std::string &msg) const;
};

#endif
