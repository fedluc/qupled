#include "util/database.hpp"
#include "schemes/qstls.hpp"
#include "util/format.hpp"

using namespace std;

namespace databaseUtil {

  void deleteBlobDataOnDisk(const string &dbName, int runId) {
    QstlsUtil::deleteBlobDataOnDisk(dbName, runId);
  }
} // namespace databaseUtil
