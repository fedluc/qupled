#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <SQLiteCpp/SQLiteCpp.h>

#include "fixtures/input_builders.hpp"
#include "fixtures/temp_paths.hpp"
#include "schemes/qstls.hpp"

TEST(QstlsApiAndUtilTest, ConstructorRejectsUnsupported2DConfiguration) {
  auto in2d = testFixtures::makeQstlsInput(
      "QSTLS", dimensionsUtil::Dimension::D2, 1.0, 0.7, 2);
  EXPECT_THROW((Qstls(in2d, false)), std::runtime_error);
}

TEST(QstlsApiAndUtilTest, ConstructorInitializesAdrFixedStorageShape) {
  auto inFinite = testFixtures::makeQstlsInput(
      "QSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.8, 1);

  Qstls solver(inFinite, false);
  const auto nx = solver.getWvg().size();
  const auto nl = static_cast<size_t>(inFinite->getNMatsubara());

  EXPECT_EQ(solver.getAdrFixed().size(0), nx);
  EXPECT_EQ(solver.getAdrFixed().size(1), nl);
  EXPECT_EQ(solver.getAdrFixed().size(2), nx);
}

TEST(QstlsApiAndUtilTest, FiniteComputePathSucceedsWithConfiguredDatabase) {

  testFixtures::ScopedTempDir tmp("qupled_qstls_compute_");
  const auto dbPath = (tmp.path() / "qstls.sqlite").string();
  {
    SQLite::Database db(dbPath, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
    db.exec("CREATE TABLE runs (id INTEGER PRIMARY KEY);");
    db.exec("INSERT INTO runs (id) VALUES (1);");
  }

  auto inFinite = testFixtures::makeQstlsInput(
      "QSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.8, 1);
  databaseUtil::DatabaseInfo dbInfo;
  dbInfo.name = dbPath;
  dbInfo.blobStorage = tmp.path().string();
  dbInfo.runId = 1;
  dbInfo.runTableName = "runs";
  inFinite->setDatabaseInfo(dbInfo);

  Qstls solver(inFinite, false);
  EXPECT_EQ(solver.compute(), 0);
  EXPECT_EQ(solver.getSsf().size(), solver.getWvg().size());
}

TEST(QstlsApiAndUtilTest, AdrGroundReturnsZeroAtXZero) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf{0.0, 0.5, 0.9, 1.0, 1.0};
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);
  QstlsUtil::AdrGround adr0(0.0, 1.0, ssfi, 2.0, itg2);
  EXPECT_DOUBLE_EQ(adr0.get(), 0.0);
}

TEST(QstlsApiAndUtilTest, SsfGroundUsesZeroCouplingFallback) {
  const std::vector<double> wvg{0.0, 0.5, 1.0, 1.5, 2.0};
  const std::vector<double> ssf{0.0, 0.5, 0.9, 1.0, 1.0};
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto in = testFixtures::makeBaseInput(
      "QSTLS", dimensionsUtil::Dimension::D3, 0.0, 0.0, 2);
  auto itg1 = std::make_shared<Integrator1D>(1.0e-8);
  QstlsUtil::SsfGround ssfGround(1.0, 0.8, 2.0, ssfi, itg1, in);
  EXPECT_NEAR(ssfGround.get(), 0.8, 1.0e-12);
}

TEST(QstlsApiAndUtilTest, DeleteBlobDataOnDiskRemovesRegisteredFiles) {
  testFixtures::ScopedTempDir tmp("qupled_qstls_blob_cleanup_");
  const auto dbPath = (tmp.path() / "test.sqlite").string();

  {
    SQLite::Database db(dbPath, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
    db.exec("CREATE TABLE fixed (run_id INTEGER NOT NULL, name TEXT NOT NULL, "
            "value TEXT NOT NULL);");

    const auto f1 = tmp.path() / "blob1.bin";
    const auto f2 = tmp.path() / "blob2.bin";
    std::ofstream(f1.string()).put('a');
    std::ofstream(f2.string()).put('b');

    SQLite::Statement s1(
        db, "INSERT INTO fixed (run_id, name, value) VALUES (?, ?, ?);");
    s1.bind(1, 42);
    s1.bind(2, "a");
    s1.bind(3, f1.string());
    s1.exec();

    SQLite::Statement s2(
        db, "INSERT INTO fixed (run_id, name, value) VALUES (?, ?, ?);");
    s2.bind(1, 42);
    s2.bind(2, "b");
    s2.bind(3, f2.string());
    s2.exec();
  }

  QstlsUtil::deleteBlobDataOnDisk(dbPath, 42);
  EXPECT_FALSE(std::filesystem::exists(tmp.path() / "blob1.bin"));
  EXPECT_FALSE(std::filesystem::exists(tmp.path() / "blob2.bin"));
}
