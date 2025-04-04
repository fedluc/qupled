#include "num_util.hpp"
#include <string>

struct DatabaseInfo {
  // Constructors
  DatabaseInfo(const std::string &name_,
               const int runId_,
               const std::string &runTableName_)
      : name(name_),
        runId(runId_),
        runTableName(runTableName_) {}
  DatabaseInfo()
      : DatabaseInfo("", numUtil::iNaN, "") {}
  // Database name
  const std::string name;
  // Run id in the database
  const int runId;
  // Name of the table with the runs in the database
  const std::string runTableName;
};