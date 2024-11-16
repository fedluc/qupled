#include "logger.hpp"

using namespace std;

void Logger::print(const string &msg) const {
  if (verbose) { cout << msg; }
}

void Logger::println(const string &msg) const {
  if (verbose) { cout << msg << endl; }
}
