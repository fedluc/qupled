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

Installing qupled
-----------------

External dependencies
~~~~~~~~~~~~~~~~~~~~~

The installation of qupled may require compiling some C++ code, depending on the platform and installation method.
Therefore, ensure the following dependencies are met before attempting to install or run qupled:

  - `CMake <https://cmake.org/download/>`_
  - `Boost <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_
  - `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
  - `OpenMP <https://en.wikipedia.org/wiki/OpenMP>`_
  - `Open-MPI <https://www.open-mpi.org/software/ompi/v5.0/>`_

For linux distributions all these dependencies can be installed with

.. code-block:: console

   sudo apt-get install -y cmake libboost-all-dev libopenmpi-dev libgsl-dev libomp-dev libfmt-dev python3-dev

For macOS they can be installed directly from homebrew

.. code-block:: console

   brew install cmake gsl libomp openmpi fmt boost-python3

Install with pip
~~~~~~~~~~~~~~~~

Qupled can be installed as a pip package by running

.. code-block:: console

   pip install qupled
		
This will also install all the python packages that are necessary for running the package.

Install from source
~~~~~~~~~~~~~~~~~~~

Qupled and all its python dependencies can also be installed from source by running

.. code-block:: console

   git clone https://github.com/fedluc/qupled.git
   cd qupled
   pip install -r dev/requirements.txt
   ./devtool build
   ./devtool test
   ./devtool install
