# Set cmake version
cmake_minimum_required (VERSION 3.30)

# Define project name
project(qupled LANGUAGES CXX)

# Set output directory for the project
set(QUPLED_OUTPUT_DIR ${CMAKE_BINARY_DIR}/qupled)

# Add subdirectories
add_subdirectory(python)
add_subdirectory(src)

# Install
install(DIRECTORY ${QUPLED_OUTPUT_DIR}
	DESTINATION ${QUPLED_INSTALL_DIR})

# Pip package
add_custom_target(pip
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/release
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/python/pyproject.toml ${CMAKE_BINARY_DIR}/release/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/python/MANIFEST.in ${CMAKE_BINARY_DIR}/release/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/qupled ${CMAKE_BINARY_DIR}/release/qupled
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR}/release python -m build
)
add_dependencies(pip qupled)

# Formatting
add_custom_target(format
    COMMAND black ${CMAKE_SOURCE_DIR}
    COMMAND find ${CMAKE_SOURCE_DIR}/src -name "*.cpp" | xargs clang-format --style=file -i
    COMMAND find ${CMAKE_SOURCE_DIR}/include -name "*.hpp" | xargs clang-format --style=file -i
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Docs
add_custom_target(docs
    COMMAND sphinx-build -b html ${CMAKE_SOURCE_DIR}/docs docs
)
