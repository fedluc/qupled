# ------------------------------------------------------------
# Shared target helpers.
#
# These wrappers centralize repeated target settings so compile/link
# behavior stays uniform across module and test targets.
# ------------------------------------------------------------

# Apply project-wide C++ standard and warning/optimization flags.
function(native_apply_common_compile_settings target)
	target_compile_features(${target} PRIVATE cxx_std_20)
	target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -O3)
endfunction()

# Link common numerical/runtime dependencies.
function(native_link_common_deps target)
	target_link_libraries(${target} PRIVATE
		OpenMP::OpenMP_CXX
		GSL::gsl
		GSL::gslcblas
	)
endfunction()

# Link normalized SQLite dependencies exported by Dependencies.cmake.
function(native_link_sqlite_deps target)
	target_link_libraries(${target} PRIVATE
		${SQLITECPP_LINK_LIB}
		${SQLITE3_LINK_LIB}
	)
endfunction()
