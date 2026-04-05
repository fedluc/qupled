Python API
==========

**qupled** is a hybrid Python/C++ package designed to simulate and analyze dielectric response 
schemes, both classical and quantum. **Python** orchestrates the workflow: setting up simulations, 
managing input/output, and storing results in an SQLite database. **C++** performs the heavy 
numerical computations, accessed via the 
`Boost.Python <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_ bindings.

To run a simulation, users configure input parameters through Python classes 
(e.g., :obj:`qupled.schemes.hf.Input`) and launch the computation using the corresponding solver 
class (e.g., :obj:`qupled.schemes.hf.HF`). Results are returned in dedicated result objects 
(e.g., :obj:`qupled.schemes.hf.Result`) and automatically written to an SQLite database for 
post-processing or analysis.

The database is implemented in standard SQLite format and can be accessed using any 
compatible external tool (e.g., `sqlite3`, DB Browser for SQLite). For convenience, 
the package provides a dedicated Python interface in the :obj:`qupled.postprocess.output.Database` 
class, which includes helper methods to query, inspect, and extract stored simulation 
data programmatically.

The sections below describe the supported schemes in more detail, including their 
respective input, solver, and result classes.


Classic schemes
---------------

HF scheme
~~~~~~~~~~

The :obj:`qupled.schemes.hf` module is used to setup and perform all the necessary calculations
for the solution of the `Hartree-Fock approximation <https://arxiv.org/abs/2411.04904>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.hf.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.hf.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.hf.Solver
    :members:

.. autoclass:: qupled.schemes.hf.Input
    :members:
    :exclude-members: to_native

.. autoclass:: qupled.schemes.hf.Result
    :members:
    :exclude-members: compute_rdf, from_native

Rpa scheme
~~~~~~~~~~

The :obj:`qupled.schemes.rpa` module is used to setup and perform all the necessary calculations
for the solution of the `Random-Phase Approximation <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.rpa.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.hf.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.rpa.Solver
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.rpa.Input	       
    :show-inheritance:
    :members:

Stls scheme
~~~~~~~~~~~		      

The :obj:`qupled.schemes.stls` module is used to setup and perform all the necessary calculations
for the solution of the `Stls scheme <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.stls.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.stls.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.stls.Solver
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.stls.Input	       
    :show-inheritance:
    :members:
      
.. autoclass:: qupled.schemes.stls.Result
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.stls.Guess
    :members:
    :exclude-members: to_native

Stls-IET schemes
~~~~~~~~~~~~~~~~	      

The :obj:`qupled.schemes.stlsiet` module is used to setup and perform all the necessary calculations
for the solution of the `Stls-IET schemes <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.stlsiet.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.stlsiet.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.stlsiet.Solver
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.stlsiet.Input	       
    :show-inheritance:
    :members:
      
.. autoclass:: qupled.schemes.stlsiet.Result
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.stlsiet.Guess
    :members:
    :exclude-members: to_native

VSStls scheme
~~~~~~~~~~~~~

The :obj:`qupled.schemes.vsstls` module is used to setup and perform all the necessary calculations
for the solution of the `VSStls scheme <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.115123>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.vsstls.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.vsstls.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.vsstls.Solver
    :show-inheritance:
    :members:
    :exclude-members: get_free_energy_integrand

.. autoclass:: qupled.schemes.vsstls.Input	       
    :show-inheritance:
    :members:
      
.. autoclass:: qupled.schemes.vsstls.Result
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.vsstls.FreeEnergyIntegrand
    :members:
    :exclude-members: to_native

Hybrid schemes
--------------

ESA scheme
~~~~~~~~~~

The :obj:`qupled.schemes.esa` module is used to setup and perform all the necessary calculations
for the solution of the `Effective Static Approximation <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.esa.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.hf.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.esa.Solver
    :show-inheritance:
    :members:

.. autoclass:: qupled.schemes.esa.Input	       
    :show-inheritance:
    :members:

Quantum schemes
---------------

Qstls scheme
~~~~~~~~~~~~

The :obj:`qupled.schemes.qstls` module is used to setup and perform all the necessary calculations
for the solution of the `Qstls scheme <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_.
TThe solution parameters are specified with a dedicated class called :obj:`qupled.schemes.qstls.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.stls.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.qstls.Solver
    :show-inheritance:
    :members:
    :exclude-members: find_fixed_adr_in_database

.. autoclass:: qupled.schemes.qstls.Input	       
    :show-inheritance:
    :members:

Qstls-IET scheme
~~~~~~~~~~~~~~~~

The :obj:`qupled.schemes.qstlsiet` module is used to setup and perform all the necessary calculations
for the solution of the `Qstls-IET schemes <https://pubs.aip.org/aip/jcp/article/158/14/141102/
2877795/Quantum-version-of-the-integral-equation-theory>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.schemes.qstlsiet.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.schemes.stlsiet.Result` 
and written to the output database.

.. autoclass:: qupled.schemes.qstlsiet.Solver
    :show-inheritance:
    :members:
    :exclude-members: find_fixed_adr_in_database

.. autoclass:: qupled.schemes.qstlsiet.Input	       
    :show-inheritance:
    :members:

QVSStls scheme
~~~~~~~~~~~~~~~~

The :obj:`qupled.schemes.qvsstls.QVSStls` class is used to setup and perform all the necessary calculations
for the solution of the QVStls schemes. The solution parameters are specified with a dedicated class 
called :obj:`qupled.schemes.qvsstls.Input`. After the solution is completed the results are stored in an 
object :obj:`qupled.schemes.vsstls.Result` and written to the output database.

.. autoclass:: qupled.schemes.qvsstls.Solver
    :show-inheritance:
    :members:
    :exclude-members: get_free_energy_integrand

.. autoclass:: qupled.schemes.qvsstls.Input	       
    :show-inheritance:
    :members:

Finite Size Correction
----------------------

The :obj:`qupled.postprocess.finite_size_correction.FiniteSizeCorrection` class is used to setup and perform 
all the necessary calculations to compute the finite size correction. The solution parameters are 
specified with a dedicated class called :obj:`qupled.postprocess.finite_size_correction.Input`. After the 
solution is completed the results are stored in an object :obj:`qupled.postprocess.finite_size_correction.Result`
and written to the output database.

.. autoclass:: qupled.postprocess.finite_size_correction.FiniteSizeCorrection
    :members:

.. autoclass:: qupled.postprocess.finite_size_correction.Input       
    :members:
      
.. autoclass:: qupled.postprocess.finite_size_correction.Result
    :members:

Output database
---------------

The output generated by qupled is stored in an SQL database named `qupled.db`, which is created in 
the directory where qupled is executed. If the database doesn't already exist, it will be created 
automatically; otherwise, new results are added to the existing database. The :obj:`qupled.postprocess.output.DataBase`
class provides functionality for accessing the data stored in the database.

.. autoclass:: qupled.postprocess.output.DataBase
    :members:

.. autoclass:: qupled.postprocess.output.OutputType
    :members:

Dimenions
---------

.. autoclass:: qupled.util.dimension.Dimension
    :members:

