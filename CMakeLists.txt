cmake_minimum_required (VERSION 3.24)

# Define the project
project(qupled LANGUAGES CXX)

# Set environment variables
set(GSL_PATH /opt/homebrew/lib)
set(OMP_PATH ${GSL_PATH})
set(ARGP_PATH ${GSL_PATH})

# Setup python environment
set(PYTHON_ROOT /opt/homebrew/opt/python@3.11/Frameworks/Python.framework/Versions/3.11)
set(PYTHON_LIBRARY ${PYTHON_ROOT}/lib/libpython3.11.dylib)
set(PYTHON_INCLUDE_DIR ${PYTHON_ROOT}/include/python3.11)
find_package(PythonInterp 3)
if (PYTHONINTERP_FOUND)
   find_package(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
   find_package(PythonLibs 3 REQUIRED)
else()
   message("Python not found")
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_CXX_STANDARD 14)

# Source code
python_add_module(qupled
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
			   ${PYTHON_INCLUDE_DIRS}
			   include)

# External libraries
find_package(OpenMP REQUIRED)
find_library(GSL_LIBRARY NAMES gsl PATHS ${GSL_PATH} REQUIRED)
find_library(GSLCBLAS_LIBRARY NAMES gslcblas PATH ${GSL_PATH} REQUIRED)
find_library(MATH_LIBRARY NAMES m REQUIRED)
find_library(ARGP_LIBRARY NAMES argp ${ARGP_PATH} REQUIRED)

# Link external libraries
target_link_libraries(qupled
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
