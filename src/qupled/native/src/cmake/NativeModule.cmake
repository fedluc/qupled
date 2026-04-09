# ------------------------------------------------------------
# Native Python extension target definition.
#
# This module assembles the pybind extension and applies common
# compile/link settings plus platform-specific link dependencies.
# ------------------------------------------------------------

# Main extension target consumed by Python.
pybind11_add_module(native ${NATIVE_SOURCES})

# Project include path required by native C++ sources.
target_include_directories(native PRIVATE
	../include
)

# Shared project defaults from TargetHelpers.
native_apply_common_compile_settings(native)
native_link_common_deps(native)
native_link_sqlite_deps(native)

# Optional MPI linkage (only when USE_MPI is enabled).
if(USE_MPI)
	target_link_libraries(native PRIVATE MPI::MPI_CXX)
endif()

# Apple-only formatting dependency.
if(APPLE)
	target_link_libraries(native PRIVATE fmt::fmt)
endif()
