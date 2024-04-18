# Set cmake version
cmake_minimum_required (VERSION 3.24)

# Define project name
project(qupled LANGUAGES CXX)

# Set output directory for the project
set(QUPLED_OUTPUT_DIR ${CMAKE_BINARY_DIR}/qupled)

find_package(Boost COMPONENTS python3 numpy3 REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_DIRS})
endif()

find_package(GSL REQUIRED)

if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
endif()

# Add subdirectories
add_subdirectory(python)
add_subdirectory(src)

# Install
install(DIRECTORY ${QUPLED_OUTPUT_DIR}
	DESTINATION ${QUPLED_INSTALL_DIR})
