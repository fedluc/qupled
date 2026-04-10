# ------------------------------------------------------------
# Source list registry.
#
# The core list is reused by both the Python extension and test-support
# library to keep behavior consistent and avoid duplication.
# ------------------------------------------------------------

# Production C++ implementation units (no Python glue)
set(NATIVE_CORE_SOURCES
	util/logger.cpp
	util/database.cpp
	util/numerics.cpp
	util/vector2D.cpp
	util/vector3D.cpp
	util/mpi_util.cpp
	util/num_util.cpp
	util/vector_util.cpp
	util/dimensions_util.cpp
	thermo/internal_energy.cpp
	thermo/free_energy.cpp
	thermo/rdf.cpp
	thermo/itcf.cpp
	thermo/thermo_util.cpp
	thermo/chemical_potential.cpp
	vs/vsmanager.cpp
	vs/vsbase.cpp
	schemes/input.cpp
	schemes/hf.cpp
	schemes/rpa.cpp
	schemes/esa.cpp
	schemes/stls.cpp
	schemes/iet.cpp
	schemes/stlsiet.cpp
	schemes/qstls.cpp
	schemes/qstlsiet.cpp
	schemes/vsstls.cpp
	schemes/qvsstls.cpp
)

# Python interface translation units.
# Separated so native tests do not depend on Python-binding sources.
set(NATIVE_PYTHON_SOURCES
	python_interface/util.cpp
	python_interface/inputs.cpp
	python_interface/schemes.cpp
	python_interface/utilities.cpp
	python_interface/native.cpp
)

# Final source list used by `pybind11_add_module(native ...)`.
set(NATIVE_SOURCES
	${NATIVE_CORE_SOURCES}
	${NATIVE_PYTHON_SOURCES}
)
