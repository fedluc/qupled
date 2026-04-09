#include <gtest/gtest.h>

#include <gsl/gsl_errno.h>

#include "util/mpi_util.hpp"

namespace {
class MpiEnvironment : public ::testing::Environment {
public:
  void SetUp() override {
    gsl_set_error_handler_off();
    MPIUtil::init();
  }
  void TearDown() override { MPIUtil::finalize(); }
};

[[maybe_unused]] const ::testing::Environment *const mpi_env =
    ::testing::AddGlobalTestEnvironment(new MpiEnvironment());
} // namespace
