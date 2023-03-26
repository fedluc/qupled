cmake_minimum_required (VERSION 3.24)

# Define the project
project(qupled LANGUAGES CXX)

# Set environment variables
set(GSL_PATH /opt/homebrew/lib)
set(OMP_PATH ${GSL_PATH})
set(ARGP_PATH ${GSL_PATH})

# Set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_CXX_STANDARD 14)

# Source code
add_executable(qupled
	       src/chemicalpotential.cpp
	       src/exec.cpp
	       src/input.cpp
	       src/numerics.cpp
	       src/qstls.cpp
	       src/stls.cpp
	       src/util.cpp)

# Set include directories
target_include_directories(qupled PUBLIC include)

# External libraries
find_package(OpenMP REQUIRED)
find_library(GSL NAMES gsl PATHS ${GSL_PATH} REQUIRED)
find_library(GSLCBLAS NAMES gslcblas PATH ${GSL_PATH} REQUIRED)
find_library(MATH NAMES m REQUIRED)
find_library(ARGP NAMES argp ${ARGP_PATH} REQUIRED)

# Link external libraries
target_link_libraries(qupled
		      OpenMP::OpenMP_CXX
		      ${GSL}
		      ${GSLCBLAS}
		      ${MATH}
		      ${OMP}
		      ${ARGP})

			   
