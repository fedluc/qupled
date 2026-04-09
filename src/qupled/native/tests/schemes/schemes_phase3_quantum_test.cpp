#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

#include <SQLiteCpp/SQLiteCpp.h>

#include "fixtures/input_builders.hpp"
#include "fixtures/temp_paths.hpp"
#include "schemes/qstls.hpp"
#include "schemes/qstlsiet.hpp"
#include "schemes/qvsstls.hpp"

TEST(SchemesPhase3Test, QuantumConstructorsValidateUnsupportedConfigurations) {
  auto qstls_2d =
      testFixtures::makeQstlsInput("QSTLS", dimensionsUtil::Dimension::D2, 1.0, 0.8, 2);
  EXPECT_THROW(Qstls(qstls_2d, false), std::runtime_error);

  auto qstlsiet_ground =
      testFixtures::makeQstlsIetInput("QSTLS-HNC", "sqrt", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  EXPECT_THROW({ QstlsIet tmp(qstlsiet_ground); }, std::runtime_error);

  auto qvs_ground =
      testFixtures::makeQVSStlsInput("QVSSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.0, 2);
  EXPECT_THROW({ QVSStls tmp(qvs_ground); }, std::runtime_error);

  auto qstls_ok =
      testFixtures::makeQstlsInput("QSTLS", dimensionsUtil::Dimension::D3, 1.0, 0.8, 2);
  EXPECT_NO_THROW(Qstls solver(qstls_ok, false));
}

TEST(SchemesPhase3Test, DeleteBlobDataOnDiskRemovesRegisteredFiles) {
  testFixtures::ScopedTempDir tmp("qupled_qstls_blob_cleanup_");
  const auto db_path = (tmp.path() / "test.sqlite").string();

  {
    SQLite::Database db(db_path, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
    db.exec("CREATE TABLE fixed (run_id INTEGER NOT NULL, name TEXT NOT NULL, value TEXT NOT NULL);");

    const auto f1 = tmp.path() / "blob1.bin";
    const auto f2 = tmp.path() / "blob2.bin";
    std::ofstream(f1.string()).put('a');
    std::ofstream(f2.string()).put('b');

    SQLite::Statement s1(db, "INSERT INTO fixed (run_id, name, value) VALUES (?, ?, ?);");
    s1.bind(1, 42);
    s1.bind(2, "adr_a");
    s1.bind(3, f1.string());
    s1.exec();

    SQLite::Statement s2(db, "INSERT INTO fixed (run_id, name, value) VALUES (?, ?, ?);");
    s2.bind(1, 42);
    s2.bind(2, "adr_b");
    s2.bind(3, f2.string());
    s2.exec();

    SQLite::Statement s3(db, "INSERT INTO fixed (run_id, name, value) VALUES (?, ?, ?);");
    s3.bind(1, 7);
    s3.bind(2, "other_run");
    s3.bind(3, (tmp.path() / "other.bin").string());
    s3.exec();
  }

  QstlsUtil::deleteBlobDataOnDisk(db_path, 42);

  EXPECT_FALSE(std::filesystem::exists(tmp.path() / "blob1.bin"));
  EXPECT_FALSE(std::filesystem::exists(tmp.path() / "blob2.bin"));
}

TEST(SchemesPhase3Test, QAdderReturnsNearZeroForUnitStructureFactor) {
  const std::vector<double> wvg{0.1, 0.4, 0.8, 1.2, 1.6};
  const std::vector<double> ssf(wvg.size(), 1.0);
  auto ssfi = std::make_shared<Interpolator1D>(wvg, ssf);
  auto itg1 = std::make_shared<Integrator1D>(1.0e-8);
  auto itg2 = std::make_shared<Integrator2D>(1.0e-6);

  const QAdder qadder(0.7, 0.0, wvg.front(), wvg.back(), {}, itg1, itg2, ssfi);
  EXPECT_NEAR(qadder.get(), 0.0, 1.0e-6);
}
