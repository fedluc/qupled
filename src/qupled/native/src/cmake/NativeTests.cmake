# ------------------------------------------------------------
# Native C++ unit test targets.
#
# This module is included only when BUILD_NATIVE_TESTS=ON and defines:
# - GoogleTest acquisition/discovery
# - a shared test-support library based on production core sources
# - the native test executable registered with CTest
# ------------------------------------------------------------

include(CTest)
enable_testing()

# Test framework setup.
# Prefer system/package-manager GTest; fetch from source as fallback.
find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
	FetchContent_Declare(
		googletest
		GIT_REPOSITORY https://github.com/google/googletest.git
		GIT_TAG v1.14.0
	)
	FetchContent_MakeAvailable(googletest)
endif()

# Gather all test translation units under ../tests.
file(GLOB_RECURSE NATIVE_TEST_FILES CONFIGURE_DEPENDS
	../tests/*.cpp
)

# Shared support library for tests.
# This keeps test binaries aligned with production core logic while
# excluding Python-binding translation units.
add_library(native_test_support STATIC ${NATIVE_CORE_SOURCES})

target_include_directories(native_test_support PRIVATE
	../include
)

native_apply_common_compile_settings(native_test_support)
native_link_common_deps(native_test_support)
native_link_sqlite_deps(native_test_support)

# Unit test executable (all discovered tests).
add_executable(native_tests
	${NATIVE_TEST_FILES}
)

target_include_directories(native_tests PRIVATE
	../include
	../tests
)

native_apply_common_compile_settings(native_tests)

target_link_libraries(native_tests PRIVATE
	GTest::gtest_main
	native_test_support
)
native_link_common_deps(native_tests)

include(GoogleTest)
# Register each GoogleTest case with CTest for `ctest` execution.
gtest_discover_tests(native_tests)
