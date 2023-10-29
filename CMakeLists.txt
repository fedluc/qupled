cmake_minimum_required (VERSION 3.24)

set(CMAKE_DEBUG_POSTFIX _debug)

# Define the project
project(qupled LANGUAGES CXX)

# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Boost
	     COMPONENTS python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED
	     COMPONENTS numpy3 REQUIRED)

# Set compiler flags
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -O3)

# Source code
Python3_add_library(qupled MODULE
		    src/chemicalpotential.cpp
		    src/input.cpp
		    src/numerics.cpp
		    src/python_modules.cpp
		    src/qstls.cpp
		    src/stls.cpp
		    src/util.cpp
		    src/vsstls.cpp)

# Include directories
target_include_directories(qupled PUBLIC
			   ${Boost_INCLUDE_DIRS}
			   ${Python3_INCLUDE_DIRS}
			   include)

# External libraries
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Link external libraries
target_link_libraries(qupled PUBLIC
		      OpenMP::OpenMP_CXX
		      ${Python3_LIBRARIES}
		      ${Boost_LIBRARIES}
		      GSL::gsl
		      GSL::gslcblas)		   

# Create python package							
add_custom_target(copy_python COMMAND
		  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/copy-python.cmake)
add_dependencies(qupled copy_python)
set_target_properties(qupled PROPERTIES LIBRARY_OUTPUT_DIRECTORY qupled)

# Install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/qupled
	DESTINATION ${Python3_SITELIB})

# Clean target
add_custom_target(true-clean COMMAND
		  COMMAND ${CMAKE_BUILD_TOOL} clean
		  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/true-clean.cmake)

