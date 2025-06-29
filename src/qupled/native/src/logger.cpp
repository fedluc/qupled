#include "logger.hpp"
#include "mpi_util.hpp"

using namespace std;

bool Logger::log_message() const { return verbose && MPIUtil::isRoot(); }

void Logger::print(const string &msg) const {
  if (log_message()) { cout << msg; }
}

void Logger::println(const string &msg) const {
  if (log_message()) { cout << msg << endl; }
}
