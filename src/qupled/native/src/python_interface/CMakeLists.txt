# -------------------------------------
# Python Interface Object Library
# -------------------------------------
add_library(native_python_interface OBJECT
    python_util.cpp
    python_wrappers.cpp
    python_modules.cpp
)

# -------------------------------------
# Target Definition
# -------------------------------------
# Include directories
target_include_directories(native_python_interface PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/python_interface
)

# Link external libraries
target_link_libraries(native_python_interface PUBLIC
    Boost::python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}
    Boost::numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}
    Python3::Python
)