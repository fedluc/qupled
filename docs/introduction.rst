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

Installation
------------

Install with pip
~~~~~~~~~~~~~~~~

Qupled can be installed as a pip package by running

.. code-block:: console

   pip install qupled

This will also install the Python packages that are necessary for running
qupled. Depending on your platform and installation method, you may also need
the :ref:`runtime_dependencies`.

Install from source
~~~~~~~~~~~~~~~~~~~

If you want full control over your qupled installation, clone the repository:

.. code-block:: console

   git clone https://github.com/fedluc/qupled.git
   cd qupled

Then install the :ref:`source_build_dependencies` for your platform. Once they
are available, install the package from the repository root:

.. code-block:: console

   curl -LsSf https://astral.sh/uv/install.sh | sh
   uvx foga install --target qupled

Install with MPI support
~~~~~~~~~~~~~~~~~~~~~~~~

MPI support is built from source. After installing the
:ref:`source_build_dependencies`, including the MPI dependencies, run from the
repository root:

.. code-block:: console

   uvx foga install --target qupled --profile mpi

Architecture
------------

Qupled is a hybrid Python/C++ package. The Python layer handles the user-facing API, input
validation, and result storage, while the C++ native library performs the heavy numerical
computations. The two layers communicate through `pybind11 <https://pybind11.readthedocs.io>`_
bindings compiled into the ``qupled.native`` extension module.

.. mermaid::

   %%{init: {'theme': 'default', 'themeVariables': {'primaryColor': '#dbeafe', 'primaryBorderColor': '#2563eb', 'primaryTextColor': '#1e3a5f', 'edgeLabelBackground': '#f0f4ff', 'lineColor': '#2563eb'}}}%%
   graph TD
       User(["User"])
       Input["Input class<br/>(qupled.schemes.*)"]
       Solver["Solver class<br/>(qupled.schemes.*)"]
       Bindings["pybind11 bindings<br/>(qupled.native)"]
       Cpp["C++ native library<br/>(schemes / thermo / util)"]
       DB[("SQLite database<br/>qupled.db")]
       Post["Post-processing<br/>(qupled.postprocess)"]

       User -->|"configures"| Input
       User -->|"calls compute()"| Solver
       Solver -->|"passes input to"| Bindings
       Bindings -->|"invokes"| Cpp
       Cpp -->|"returns results"| Bindings
       Bindings -->|"returns"| Solver
       Solver -->|"writes results to"| DB
       DB -->|"read via"| Post

Units
-----

All the calculations are performed in normalized units. The wave vectors are normalized to the
Fermi wave-vector and the frequencies are normalized to :math:`2\pi E_{\mathrm{f}}/h`. Here :math:`E_{\mathrm{f}}`
is the Fermi energy and :math:`h` is Planck's constant.

Limitations
-----------

* Quantum schemes are currently only available for 3D systems; 2D calculations
  are not implemented for QSTLS, QSTLS-IET, and QVS.
* Ground state (zero temperature) calculations are not available for the
  QSTLS-IET and QVS schemes.
