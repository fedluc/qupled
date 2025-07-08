#include "database.hpp"
#include "format.hpp"
#include "qstls.hpp"
#include <SQLiteCpp/SQLiteCpp.h>

using namespace std;

namespace databaseUtil {

  void deleteBlobDataOnDisk(const DatabaseInfo &dbInfo, int runId) {
    try {
      // Setup the database connection
      SQLite::Database db(dbInfo.name, SQLite::OPEN_READONLY);
      const string select = formatUtil::format(QstlsUtil::SQL_SELECT_RUN_ID,
                                               QstlsUtil::SQL_TABLE_NAME);
      SQLite::Statement statement(db, select);
      statement.bind(1, runId);
      std::vector<std::string> filePaths;
      // Execute and collect file paths
      while (statement.executeStep()) {
        filePaths.push_back(statement.getColumn(0).getString());
      }
      // Delete each file path
      for (const auto &path : filePaths) {
        try {
          if (std::filesystem::remove(path)) {
            // Successfully deleted
          } else {
            std::cerr << "Warning: file not found: " << path << std::endl;
          }
        } catch (const std::filesystem::filesystem_error &e) {
          std::cerr << "Warning: failed to delete " << path << ": " << e.what()
                    << std::endl;
        }
      }
    } catch (const SQLite::Exception &e) {
      std::cerr << "SQLite error: " << e.what() << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "General error: " << e.what() << std::endl;
    }
  }
} // namespace databaseUtil