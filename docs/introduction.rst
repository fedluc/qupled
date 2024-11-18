Overview
========

Qupled is a package that can be used to compute the properties of quantum one component
plasmas via theoretical approaches based on the dielectric formalism. The theoretical
approaches which can be solved with qupled include:

  * The classical RPA scheme (`Böhm and Pines <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_)
  * The classical STLS scheme (`Tanaka and Ichimaru <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_)
  * The classical STLS-HNC scheme (`Tanaka <https://pubs.aip.org/aip/jcp/article/145/21/214104/196066/Correlational-and-thermodynamic-properties-of>`_)
  * The classical STLS-IET scheme (`Tolias and collaborators <https://pubs.aip.org/aip/jcp/article/155/13/134115/353165/Integral-equation-theory-based-dielectric-scheme>`_)
  * The classical VS-STLS scheme (`Vashishta and Singwi <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875>`_)
  * The hybrid ESA scheme (`Dornheim and collaborators <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.165102>`_)
  * The quantum STLS (QSTLS) scheme (`Schweng and Böhm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_)
  * The quantum QSTLS-IET scheme (`Tolias and others <https://pubs.aip.org/aip/jcp/article/158/14/141102/2877795/Quantum-version-of-the-integral-equation-theory>`_)
  * The quantum VS-STLS (QVS) scheme

Qupled supports both MPI and OpenMP parallelizations to handle the most computationally-intensive
calculations in the quantum and in the classical VS-STLS scheme.
    
Limitations
-----------

Ground state (zero temperature) calculations are not available for the quantum schemes (QSTLS, QSTLS-IET and QVS).

Units
-----

All the calculations are performed in normalized units. The wave vectors are normalized to the
Fermi wave-vector and the frequencies are normalized to :math:`2\pi E_{\mathrm{f}}/h`. Here :math:`E_{\mathrm{f}}`
is the Fermi energy and :math:`h` is Planck's constant.

Installing qupled
-----------------

Qupled can be installed as a pip package by running

.. code-block:: console

   pip install qupled
		
This will also install all the python packages that are necessary for running the package. Note that, depending on the platform,
the installation procedure might rquire to compile some C++ code. Hence, before trying to install or run qupled make sure that
the following dependencies are satisfied:

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
