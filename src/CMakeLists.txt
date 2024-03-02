# Setup python environment
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Boost
	     COMPONENTS python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED
	     COMPONENTS numpy3 REQUIRED)

# External libraries
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Set compiler flags
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -O3)

# Source code
Python3_add_library(qupled MODULE
		    chemical_potential.cpp
		    input.cpp
		    numerics.cpp
		    qstls.cpp
		    rpa.cpp
		    stls.cpp
		    util.cpp
		    vsstls.cpp
		    esa.cpp
		    qvs.cpp
		    python_wrappers.cpp
		    python_modules.cpp)

# Include directories
target_include_directories(qupled PUBLIC
			   ${Boost_INCLUDE_DIRS}
			   ${Python3_INCLUDE_DIRS}
			   ${MPI_CXX_INLUDE_PATH}
			   ../include)

# Link external libraries
target_link_libraries(qupled PUBLIC
		      OpenMP::OpenMP_CXX
		      ${Python3_LIBRARIES}
		      ${Boost_LIBRARIES}
		      ${MPI_CXX_LIBRARIES}
		      GSL::gsl
		      GSL::gslcblas)		   

# Create python package							
set_target_properties(qupled PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${QUPLED_OUTPUT_DIR})

# Define installation directory
set(QUPLED_INSTALL_DIR ${Python3_SITELIB}  CACHE INTERNAL "")
