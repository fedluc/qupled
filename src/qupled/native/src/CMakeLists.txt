# Set cmake version
cmake_minimum_required(VERSION 3.28)

# Define project name
project(native LANGUAGES CXX)

# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Development.Module)

# Setup pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# Setup fmt (only for macOS)
if(APPLE)
  message(STATUS "Configuring fmt for macOS")
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1
  )
  FetchContent_MakeAvailable(fmt)
endif()

# Setup MPI environment
option(USE_MPI "Compile with MPI support" ON)
if(USE_MPI)
	set(MPI_CXX_SKIP_MPICXX TRUE)
	find_package(MPI REQUIRED)
	add_compile_definitions(USE_MPI)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# Additional external libraries
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(SQLiteCpp REQUIRED)

# Set compiler flags
set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Wextra -Wpedantic -O3)

# Source code
pybind11_add_module(native
	logger.cpp
	database.cpp
	numerics.cpp
	vector2D.cpp
	vector3D.cpp
	internal_energy.cpp
	free_energy.cpp
	rdf.cpp
	mpi_util.cpp
	num_util.cpp
	vector_util.cpp
	thermo_util.cpp
	dimensions_util.cpp
	input.cpp
	chemical_potential.cpp
	hf.cpp
	rpa.cpp
	esa.cpp
	stls.cpp
	iet.cpp
	stlsiet.cpp
	qstls.cpp
	qstlsiet.cpp
	vsbase.cpp
	vsstls.cpp
	qvsstls.cpp
	python_interface/util.cpp
	python_interface/inputs.cpp
	python_interface/schemes.cpp
	python_interface/utilities.cpp
	python_interface/native.cpp
)

# Include directories
target_include_directories(native PUBLIC
	../include
)

# Link external libraries
target_link_libraries(native PUBLIC
	OpenMP::OpenMP_CXX
	GSL::gsl
	GSL::gslcblas
	SQLiteCpp
	SQLite::SQLite3
)
if(USE_MPI)
	target_link_libraries(native PUBLIC MPI::MPI_CXX)
endif()
if(APPLE)
	target_link_libraries(native PUBLIC fmt::fmt)
endif()

# Install
install(TARGETS native
	DESTINATION .
)