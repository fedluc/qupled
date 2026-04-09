# Set cmake version
cmake_minimum_required(VERSION 3.28)

# Define project name
project(native LANGUAGES CXX)

# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Development.Module)

# Setup pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# Setup fmt (only for macOS)
if(APPLE)
  message(STATUS "Configuring fmt for macOS")
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1
  )
  FetchContent_MakeAvailable(fmt)
endif()

# Setup MPI environment
option(USE_MPI "Compile with MPI support" ON)
option(BUILD_NATIVE_TESTS "Build native C++ unit tests with GoogleTest" OFF)
if(USE_MPI)
	set(MPI_CXX_SKIP_MPICXX TRUE)
	find_package(MPI REQUIRED)
	add_compile_definitions(USE_MPI)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# Additional external libraries
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(SQLiteCpp REQUIRED)

# Resolve SQLite targets across different CMake/package variants.
if(TARGET SQLite::SQLite3)
	set(SQLITE3_LINK_LIB SQLite::SQLite3)
elseif(TARGET SQLite3::SQLite3)
	set(SQLITE3_LINK_LIB SQLite3::SQLite3)
elseif(SQLite3_LIBRARIES)
	set(SQLITE3_LINK_LIB ${SQLite3_LIBRARIES})
else()
	message(FATAL_ERROR "Could not resolve a SQLite3 link target/library")
endif()

if(TARGET SQLiteCpp::SQLiteCpp)
	set(SQLITECPP_LINK_LIB SQLiteCpp::SQLiteCpp)
elseif(TARGET SQLiteCpp)
	set(SQLITECPP_LINK_LIB SQLiteCpp)
elseif(SQLiteCpp_LIBRARIES)
	set(SQLITECPP_LINK_LIB ${SQLiteCpp_LIBRARIES})
else()
	message(FATAL_ERROR "Could not resolve a SQLiteCpp link target/library")
endif()

# Source code
set(NATIVE_SOURCES
	util/logger.cpp
	util/database.cpp
	util/numerics.cpp
	util/vector2D.cpp
	util/vector3D.cpp
	util/mpi_util.cpp
	util/num_util.cpp
	util/vector_util.cpp
	util/dimensions_util.cpp
	thermo/internal_energy.cpp
	thermo/free_energy.cpp
	thermo/rdf.cpp
	thermo/itcf.cpp
	thermo/thermo_util.cpp
	thermo/chemical_potential.cpp
	vs/vsmanager.cpp
	vs/vsbase.cpp
	schemes/input.cpp
	schemes/hf.cpp
	schemes/rpa.cpp
	schemes/esa.cpp
	schemes/stls.cpp
	schemes/iet.cpp
	schemes/stlsiet.cpp
	schemes/qstls.cpp
	schemes/qstlsiet.cpp
	schemes/vsstls.cpp
	schemes/qvsstls.cpp

	python_interface/util.cpp
	python_interface/inputs.cpp
	python_interface/schemes.cpp
	python_interface/utilities.cpp
	python_interface/native.cpp
)

pybind11_add_module(native ${NATIVE_SOURCES})

# Include directories
target_include_directories(native PRIVATE
	../include
)

target_compile_features(native PRIVATE cxx_std_20)
target_compile_options(native PRIVATE -Wall -Wextra -Wpedantic -O3)

# Link external libraries
target_link_libraries(native PRIVATE
	OpenMP::OpenMP_CXX
	GSL::gsl
	GSL::gslcblas
	${SQLITECPP_LINK_LIB}
	${SQLITE3_LINK_LIB}
)
if(USE_MPI)
	target_link_libraries(native PRIVATE MPI::MPI_CXX)
endif()
if(APPLE)
	target_link_libraries(native PRIVATE fmt::fmt)
endif()

# Native C++ unit tests
if(BUILD_NATIVE_TESTS)
	include(CTest)
	enable_testing()

	find_package(GTest CONFIG QUIET)
	if(NOT GTest_FOUND)
		FetchContent_Declare(
		  googletest
		  GIT_REPOSITORY https://github.com/google/googletest.git
		  GIT_TAG v1.14.0
		)
		FetchContent_MakeAvailable(googletest)
	endif()

	file(GLOB_RECURSE NATIVE_TEST_FILES CONFIGURE_DEPENDS
		../tests/*.cpp
	)

	add_library(native_test_support STATIC
		util/logger.cpp
		util/numerics.cpp
		util/vector2D.cpp
		util/vector3D.cpp
		util/mpi_util.cpp
		util/num_util.cpp
		util/vector_util.cpp
		util/dimensions_util.cpp
		schemes/input.cpp
		schemes/hf.cpp
		schemes/rpa.cpp
		thermo/free_energy.cpp
		thermo/internal_energy.cpp
		thermo/rdf.cpp
		thermo/itcf.cpp
		thermo/thermo_util.cpp
		thermo/chemical_potential.cpp
		util/database.cpp
		vs/vsmanager.cpp
		vs/vsbase.cpp
		schemes/iet.cpp
		schemes/stls.cpp
		schemes/stlsiet.cpp
		schemes/esa.cpp
		schemes/qstls.cpp
		schemes/qstlsiet.cpp
		schemes/vsstls.cpp
		schemes/qvsstls.cpp
	)

	target_include_directories(native_test_support PRIVATE
		../include
	)

	target_compile_features(native_test_support PRIVATE cxx_std_20)
	target_compile_options(native_test_support PRIVATE -Wall -Wextra -Wpedantic -O3)

	target_link_libraries(native_test_support PRIVATE
		OpenMP::OpenMP_CXX
		GSL::gsl
		GSL::gslcblas
		${SQLITECPP_LINK_LIB}
		${SQLITE3_LINK_LIB}
	)

	add_executable(native_tests
		${NATIVE_TEST_FILES}
	)

	target_include_directories(native_tests PRIVATE
		../include
		../tests
	)

	target_compile_features(native_tests PRIVATE cxx_std_20)
	target_compile_options(native_tests PRIVATE -Wall -Wextra -Wpedantic -O3)

	target_link_libraries(native_tests PRIVATE
		GTest::gtest_main
		native_test_support
		OpenMP::OpenMP_CXX
		GSL::gsl
		GSL::gslcblas
	)

	include(GoogleTest)
	gtest_discover_tests(native_tests)
endif()

# Install
install(TARGETS native
	DESTINATION .
)
