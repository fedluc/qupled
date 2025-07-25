#ifndef DATABASE_HPP
#define DATABASE_HPP

#include "num_util.hpp"
#include <string>

namespace databaseUtil {

  struct DatabaseInfo {
    DatabaseInfo()
        : runId(numUtil::iNaN) {}
    // Database name
    std::string name;
    // Fodler used to store the blob data
    std::string blobStorage;
    // Run id in the database
    int runId;
    // Name of the table with the runs in the database
    std::string runTableName;
  };

  void deleteBlobDataOnDisk(const std::string &dbInfo, int runId);
  bool blobDataTableExists(const std::string &dbName);

} // namespace databaseUtil

#endif // DATABASE_HPP