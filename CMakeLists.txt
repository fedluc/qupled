# Set cmake version
cmake_minimum_required (VERSION 3.24)

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

