Overview
========

Qupled is a package that can be used to compute the properties of quantum one component
plasmas via theoretical approaches based on the dielectric formalism. The theoretical
approaches which can be solved with qupled include:

  * The classical schemes
    
    * `RPA <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_
    * `STLS <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_
    * `STLS-IET <https://pubs.aip.org/aip/jcp/article/155/13/134115/353165/Integral-equation-theory-based-dielectric-scheme>`_
    * `VS-STLS <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875>`_      
  * The quantum schemes
    
    * `QSTLS <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_
    * `QSTLS-IET <https://pubs.aip.org/aip/jcp/article/158/14/141102/2877795/Quantum-version-of-the-integral-equation-theory>`_
    * QVS
      
  * The hybrid `ESA <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.165102>`_ scheme

Qupled supports both MPI and OpenMP parallelizations to handle the most computationally-intensive
calculations in the quantum and in the classical VS-STLS scheme.
    
Limitations
-----------

Ground state (zero temperature) calculations are not available for the QSTLS-IET and QVS schemes.

Units
-----

All the calculations are performed in normalized units. The wave vectors are normalized to the
Fermi wave-vector and the frequencies are normalized to :math:`2\pi E_{\mathrm{f}}/h`. Here :math:`E_{\mathrm{f}}`
is the Fermi energy and :math:`h` is Planck's constant.

Installation
------------

Install with pip
~~~~~~~~~~~~~~~~

Qupled can be installed as a pip package by running

.. code-block:: console

   pip install qupled
		
This will also install all the python packages that are necessary for running the package. 

Install with MPI support
~~~~~~~~~~~~~~~~~~~~~~~~

If you need to use qupled with MPI support, first install the  :ref:`external_dependencies` and then run:

.. code-block:: console

   USE_MPI=ON pip install --no-binary=:all: qupled


Install from source
~~~~~~~~~~~~~~~~~~~

If you want full control over your qupled installation, you can install it also directly from the source.
Install the :ref:`external_dependencies` first, then clone the repository and run the following commands:

.. code-block:: console

   git clone https://github.com/fedluc/qupled.git
   cd qupled
   ./devtool install-deps
   ./devtool build
   ./devtool test
   ./devtool install

.. _external_dependencies:

External dependencies
~~~~~~~~~~~~~~~~~~~~~

Installing qupled may require compiling some C++ code, depending on the platform and installation method.
The following dependencies must be met before attempting to build the C++ part of qupled:

  - `CMake <https://cmake.org/download/>`_
  - `fmt <https://github.com/fmtlib/fmt>`_
  - `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
  - `OpenMP <https://en.wikipedia.org/wiki/OpenMP>`_
  - `Open-MPI <https://www.open-mpi.org/software/ompi/v5.0/>`_

**Ubuntu or Debian-based systems**

.. code-block:: console

   sudo apt-get install -y cmake libopenmpi-dev libgsl-dev libomp-dev python3-dev libsqlite3-dev libsqlitecpp-dev

**Fedora or Red Hat-based system**

.. code-block:: console

   sudo dnf install -y cmake openmpi openmpi-devel gsl-devel sqlite-devel
   cd /tmp
   git clone https://github.com/SRombauts/SQLiteCpp.git
   cd SQLiteCpp
   mkdir build && cd build 
   cmake ..
   make -j$(nproc)
   make install
   ldconfig
   cd / && rm -rf /tmp/SQLiteCpp

**macOS**

.. code-block:: console

   brew install cmake gsl libomp openmpi sqlite sqlitecpp
