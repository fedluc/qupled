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

Building qupled
---------------

Qupled is a hybrid C++/python code that can be built by using the CMake build system. The build process is
summarized as follows:

.. code-block:: console

   git clone git@github.com:fedluc/qupled.git
   cd qupled
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   cmake --build .
   pytest tests
   cmake --install .
		
This will create the python package that can be used to solve the dielectric schemes, test it and install it
in a folder accessible at runtime.  The build directory can be cleaned by running ``make true-clean``. For debugging
purposes it is also possible to build a debug configuration by typing ``cmake -DCMAKE_BUILD_TYPE=Debug ..``
followed by the build command. By default qupled is built with support for MPI, but it is also possible to compile the
code without MPI by typing ``cmake -DUSE_MPI=OFF ..`` followed by the build command.

Dependencies
------------

To build the code, the following dependencies must be satisfied:


  - `Python (v3.12) <https://www.python.org/downloads/>`_
  - `CMake (v3.26) <https://cmake.org/download/>`_
  - `Boost (v1.84) <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_
  - `GNU Scientific Library (v2.7) <https://www.gnu.org/software/gsl/>`_
  - `OpenMP (v5.0) <https://en.wikipedia.org/wiki/OpenMP>`_
  - `Open-MPI (v5.0) <https://www.open-mpi.org/software/ompi/v5.0/>`_
  - `FMT (v10.2) <https://github.com/fmtlib/fmt/>`_

To test and run the code, the following python packages must be installed

  - `matplotlib (v>=3.7) <https://matplotlib.org>`_
  - `numpy (v>=1.24)  <https://numpy.org>`_
  - `pandas (v>=2.0) <https://pandas.pydata.org>`_
  - `pytest (v>=8.0) <https://docs.pytest.org/en/8.0.x/>`_
  - `pytest-mock (v>=3.12) <https://pypi.org/project/pytest-mock/>`_
