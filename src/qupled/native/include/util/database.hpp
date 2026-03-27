#ifndef DATABASE_HPP
#define DATABASE_HPP

#include "util/num_util.hpp"
#include <string>

/**
 * @brief Utilities for SQLite database interaction used by the solvers.
 *
 * Provides a metadata struct for identifying a run's database entry and
 * helper functions for managing blob data stored alongside the database.
 */
namespace databaseUtil {

  /**
   * @brief Metadata identifying a simulation run in the SQLite database.
   *
   * Each solver run that stores results in the database is assigned a unique
   * integer @p runId within the table named @p runTableName.  Large binary
   * data (e.g., the fixed ADR component) may be stored outside the database
   * in a directory given by @p blobStorage.
   */
  struct DatabaseInfo {
    /** @brief Construct with @p runId set to the sentinel integer NaN. */
    DatabaseInfo()
        : runId(numUtil::iNaN) {}

    /** @brief Path to the SQLite database file. */
    std::string name;

    /** @brief Directory used to store large binary (blob) data on disk. */
    std::string blobStorage;

    /** @brief Integer run identifier within @p runTableName. */
    int runId;

    /** @brief Name of the table that stores per-run metadata. */
    std::string runTableName;
  };

  /**
   * @brief Delete all blob data on disk associated with a given run.
   * @param dbInfo  Path to the SQLite database file.
   * @param runId   Run identifier whose blobs should be deleted.
   */
  void deleteBlobDataOnDisk(const std::string &dbInfo, int runId);

  /**
   * @brief Check whether the blob-data table exists in a database.
   * @param dbName Path to the SQLite database file.
   * @return True if the table exists, false otherwise.
   */
  bool blobDataTableExists(const std::string &dbName);

} // namespace databaseUtil

#endif // DATABASE_HPP
