External Dependencies
=====================

Qupled uses a compiled native extension. Some external libraries are needed to
run that extension, while headers and build tools are only needed when building
qupled from source.

.. _runtime_dependencies:

Runtime Dependencies
--------------------

These libraries must be available when importing and running ``qupled.native``:

  - `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
  - `OpenMP <https://www.openmp.org/>`_
  - `SQLite <https://www.sqlite.org/>`_
  - `SQLiteCpp <https://github.com/SRombauts/SQLiteCpp>`_
  - `Open MPI <https://www.open-mpi.org/>`_ (only for MPI-enabled installs)

Package names for these runtime libraries vary by operating system and
release. If you need to install them manually, use your system package manager
to install the packages that provide the libraries listed above.

Common package-manager commands look like:

.. code-block:: console

   sudo apt-get update
   sudo apt-get install -y libgsl-dev libomp-dev libsqlite3-dev libsqlitecpp-dev
   sudo apt-get install -y openmpi-bin libopenmpi-dev  # only for MPI-enabled installs

.. code-block:: console

   brew install gsl libomp sqlite sqlitecpp
   brew install open-mpi  # only for MPI-enabled installs

.. _source_build_dependencies:

Source-Build Dependencies
-------------------------

Building qupled from source requires the runtime dependencies plus the headers
and build tools used by the native extension:

  - `CMake <https://cmake.org/download/>`_
  - `Python development headers <https://docs.python.org/3/c-api/>`_
  - `GSL development headers and libraries <https://www.gnu.org/software/gsl/>`_
  - `OpenMP development headers and libraries <https://www.openmp.org/>`_
  - `SQLite development headers and libraries <https://www.sqlite.org/>`_
  - `SQLiteCpp <https://github.com/SRombauts/SQLiteCpp>`_
  - `Open MPI development headers and libraries <https://www.open-mpi.org/>`_
    (only for MPI-enabled builds)

Preferred Installation
^^^^^^^^^^^^^^^^^^^^^^

When installing from a cloned repository, prefer the ``foga`` system dependency
target for your platform. Run these commands from the repository root:

.. code-block:: console

   curl -LsSf https://astral.sh/uv/install.sh | sh
   uvx foga install --target apt-get   # Ubuntu or Debian-based systems
   uvx foga install --target brew      # macOS

These targets intentionally install a project-oriented superset of dependencies,
including tools used for development, documentation, and tests.

Manual Installation
^^^^^^^^^^^^^^^^^^^

If ``foga`` does not provide a target for your package manager, install the
source-build dependencies manually. For example, on Fedora you need something
like:

.. code-block:: console

   sudo dnf install -y cmake gsl-devel libgomp python3-devel sqlite-devel
   cd /tmp
   git clone https://github.com/SRombauts/SQLiteCpp.git
   cd SQLiteCpp
   mkdir build && cd build
   cmake ..
   make -j$(nproc)
   sudo make install
   sudo ldconfig

Add ``openmpi`` and ``openmpi-devel`` to build with MPI support.
