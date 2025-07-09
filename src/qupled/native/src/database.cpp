#include "database.hpp"
#include "format.hpp"
#include "qstls.hpp"

using namespace std;

namespace databaseUtil {

  void deleteBlobDataOnDisk(const string &dbName, int runId) {
    QstlsUtil::deleteBlobDataOnDisk(dbName, runId);
  }
} // namespace databaseUtil
