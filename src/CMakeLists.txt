# Naming for debug target
set(CMAKE_DEBUG_POSTFIX _debug)

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
		    chemicalpotential.cpp
		    input.cpp
		    numerics.cpp
		    python_modules.cpp
		    qstls.cpp
		    rpa.cpp
		    stls.cpp
		    util.cpp
		    vsstls.cpp)

# Include directories
target_include_directories(qupled PUBLIC
			   ${Boost_INCLUDE_DIRS}
			   ${Python3_INCLUDE_DIRS}
			   ../include)

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
set_target_properties(qupled PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${QUPLED_OUTPUT_DIR})

# Define installation directory
set(QUPLED_INSTALL_DIR ${Python3_SITELIB}  CACHE INTERNAL "")
