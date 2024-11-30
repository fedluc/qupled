# Set cmake version
cmake_minimum_required (VERSION 3.30)

# Define project name
project(native LANGUAGES CXX)

# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Boost
	     COMPONENTS python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED
	     COMPONENTS numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED)


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
find_package(fmt REQUIRED)


# Set compiler flags
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -O3)


# Source code
Python3_add_library(native MODULE
		    logger.cpp
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
		    input.cpp
		    chemical_potential.cpp
		    rpa.cpp
		    esa.cpp
		    stls.cpp
		    qstls.cpp
		    vsbase.cpp
		    vsstls.cpp
		    qvs.cpp
		    python_util.cpp
		    python_wrappers.cpp
		    python_modules.cpp)


# Include directories
target_include_directories(native PUBLIC
			   ../include)


# Link external libraries
if(USE_MPI)
    target_link_libraries(native PUBLIC
			  MPI::MPI_CXX)
endif()
target_link_libraries(native PUBLIC
		      OpenMP::OpenMP_CXX
		      Boost::python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}
		      Boost::numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}
		      GSL::gsl
		      GSL::gslcblas
		      fmt::fmt)

# Install
install(TARGETS native
	DESTINATION .)
