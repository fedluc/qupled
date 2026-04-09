# ------------------------------------------------------------
# Third-party dependency bootstrap.
#
# Responsibilities of this module:
# - Configure package acquisition (find_package / FetchContent)
# - Apply platform/toolchain-specific dependency choices
# - Normalize exported target names into stable link variables
# ------------------------------------------------------------

# Python build dependency and FetchContent support
find_package(Python3 REQUIRED COMPONENTS Development.Module)
include(FetchContent)

# pybind11 is consumed by the native Python extension target.
FetchContent_Declare(
	pybind11
	GIT_REPOSITORY https://github.com/pybind/pybind11.git
	GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# fmt is required only for macOS builds in this project.
if(APPLE)
	message(STATUS "Configuring fmt for macOS")
	FetchContent_Declare(
		fmt
		GIT_REPOSITORY https://github.com/fmtlib/fmt.git
		GIT_TAG 10.2.1
	)
	FetchContent_MakeAvailable(fmt)
endif()

# Optional MPI toolchain setup.
#
# Notes:
# - `MPI_CXX_SKIP_MPICXX` avoids deprecated MPI C++ bindings.
# - The compiler is switched to the MPI wrapper when enabled.
# - A global `USE_MPI` define is added for source-level conditionals.
if(USE_MPI)
	set(MPI_CXX_SKIP_MPICXX TRUE)
	find_package(MPI REQUIRED)
	add_compile_definitions(USE_MPI)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# Core dependencies shared by production and test-support targets.
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(SQLiteCpp REQUIRED)

# Normalize SQLite target names across package variants.
# Result: SQLITE3_LINK_LIB holds a valid link item on all platforms.
if(TARGET SQLite::SQLite3)
	set(SQLITE3_LINK_LIB SQLite::SQLite3)
elseif(TARGET SQLite3::SQLite3)
	set(SQLITE3_LINK_LIB SQLite3::SQLite3)
elseif(SQLite3_LIBRARIES)
	set(SQLITE3_LINK_LIB ${SQLite3_LIBRARIES})
else()
	message(FATAL_ERROR "Could not resolve a SQLite3 link target/library")
endif()

# Normalize SQLiteCpp target names across package variants.
# Result: SQLITECPP_LINK_LIB holds a valid link item on all platforms.
if(TARGET SQLiteCpp::SQLiteCpp)
	set(SQLITECPP_LINK_LIB SQLiteCpp::SQLiteCpp)
elseif(TARGET SQLiteCpp)
	set(SQLITECPP_LINK_LIB SQLiteCpp)
elseif(SQLiteCpp_LIBRARIES)
	set(SQLITECPP_LINK_LIB ${SQLiteCpp_LIBRARIES})
else()
	message(FATAL_ERROR "Could not resolve a SQLiteCpp link target/library")
endif()
