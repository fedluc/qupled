# ------------------------------------------------------------
# Top-level entry point for the native build.
#
# This file intentionally stays small and delegates implementation
# details to module files under `cmake/`.
# ------------------------------------------------------------

# CMake and project metadata
cmake_minimum_required(VERSION 3.28)

# Define project name
project(native LANGUAGES CXX)

# Build feature toggles exposed to users/configure presets
option(USE_MPI "Compile with MPI support" ON)
option(BUILD_NATIVE_TESTS "Build native C++ unit tests with GoogleTest" OFF)

# ------------------------------------------------------------
# Core build modules (load order matters):
# 1) dependencies/toolchain knobs
# 2) source lists
# 3) shared target helper functions
# 4) native Python extension target
# ------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Dependencies.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Sources.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/TargetHelpers.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NativeModule.cmake)

# Optional native C++ test target graph
if(BUILD_NATIVE_TESTS)
	include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NativeTests.cmake)
endif()

# Install
install(TARGETS native
	DESTINATION .
)
