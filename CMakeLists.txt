cmake_minimum_required (VERSION 3.24)

# Define the project
project(qupled LANGUAGES CXX)

# Set environment variables
set(GSL_PATH /opt/homebrew/lib)
set(OMP_PATH ${GSL_PATH})
set(ARGP_PATH ${GSL_PATH})

# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Boost
	     COMPONENTS python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED
	     COMPONENTS numpy3 REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_CXX_STANDARD 14)

# Source code
Python3_add_library(qupled MODULE
		    src/chemicalpotential.cpp
		    src/input.cpp
		    src/numerics.cpp
		    src/python_modules.cpp
		    src/qstls.cpp
		    src/stls.cpp
		    src/util.cpp)

# Include directories
target_include_directories(qupled PUBLIC
			   ${Boost_INCLUDE_DIRS}
			   ${Python3_INCLUDE_DIRS}
			   include)

# External libraries
find_package(OpenMP REQUIRED)
find_library(GSL_LIBRARY NAMES gsl PATHS ${GSL_PATH} REQUIRED)
find_library(GSLCBLAS_LIBRARY NAMES gslcblas PATH ${GSL_PATH} REQUIRED)
find_library(MATH_LIBRARY NAMES m REQUIRED)
find_library(ARGP_LIBRARY NAMES argp ${ARGP_PATH} REQUIRED)

# Link external libraries
target_link_libraries(qupled PUBLIC
		      OpenMP::OpenMP_CXX
		      ${PYTHON_LIBRARIES}
		      ${Boost_LIBRARIES}
		      ${GSL_LIBRARY}
		      ${GSLCBLAS_LIBRARY}
		      ${MATH_LIBRARY}
		      ${ARGP_LIBRARY})		   

# Create python package							
add_custom_target(copy_python COMMAND
		  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/copy-python.cmake)
add_dependencies(qupled copy_python)
set_target_properties(qupled PROPERTIES LIBRARY_OUTPUT_DIRECTORY qupled)

# Clean target
add_custom_target(true-clean COMMAND
		  COMMAND ${CMAKE_BUILD_TOOL} clean
		  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/true-clean.cmake)

