The package
===========

The package is written in both Python and C++. Python can be
used to setup, run, save and postprocess a dielectric scheme of choice. C++
is used to perform the calculations involved in the solution of the schemes.
Communication between Python and C++ is handled via
`Boost <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_.
The following sections illustrate the classes used to solve classical and
quantum schemes and the classes used to handle the inputs, initial guesses and outputs.

Classic schemes
---------------

Rpa scheme
~~~~~~~~~~

The :obj:`qupled.classic.Rpa` class is used to setup and perform all the necessary calculations
for the solution of the `Random-Phase Approximation <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.Rpa.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <rpa_pandas_table>`.

.. autoclass:: qupled.classic.Rpa
    :members:
    :inherited-members:
    :exclude-members: Input, rdf


.. autoclass:: qupled.classic.Rpa.Input
    :members:

.. _rpa_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <rpa_info_table>`
     - Pandas DataFrame
     - Information on the input parameters
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function

\* available only if the radial distribution function was computed.

.. _rpa_info_table:
.. list-table:: Content of the info Dataframe
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - coupling
     - float
     - Coupling parameter
   * - degeneracy
     - float
     - Degeneracy parameter
   * - theory
     - str
     - The theory that is being solved
   * - resolution
     - float
     - Resolution in the wave-vector grid
   * - cutoff
     - float
     - Cutoff in the wave-vector grid
   * - matsubara
     - int
     - Number of matsubara frequencies

       
Stls scheme
~~~~~~~~~~~		      

The :obj:`qupled.classic.Stls` class is used to setup and perform all the necessary calculations
for the solution of the `Stls scheme <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.Stls.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <stls_pandas_table>`.

.. autoclass:: qupled.classic.Stls
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.classic.Stls.Input	       
    :members:
    :inherited-members:
    :exclude-members: print, isEqual       
       
.. _stls_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function

\* available only if the radial distribution function was computed.

.. _stls_info_table:
.. list-table:: Content of the info Dataframe
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - coupling
     - float
     - Coupling parameter
   * - degeneracy
     - float
     - Degeneracy parameter
   * - error
     - float
     - The residual error at the end of the solution
   * - theory
     - str
     - The theory that is being solved
   * - resolution
     - float
     - Resolution in the wave-vector grid
   * - cutoff
     - float
     - Cutoff in the wave-vector grid
   * - matsubara
     - int
     - Number of matsubara frequencies

Stls-IET schemes
~~~~~~~~~~~~~~~~	      

The :obj:`qupled.classic.StlsIet` class is used to setup and perform all the necessary calculations
for the solution of the `Stls-IET schemes <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.StlsIet.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <stlsiet_pandas_table>`.
    
.. autoclass:: qupled.classic.StlsIet
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.classic.StlsIet.Input 
    :members:
    :inherited-members:
    :exclude-members: print, isEqual

.. _stlsiet_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - bf
     - ndarray
     - The bridge function adder
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function

VSStls scheme
~~~~~~~~~~~~~

The :obj:`qupled.classic.VSStls` class is used to setup and perform all the necessary calculations
for the solution of the `VSStls scheme <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.115123>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.VSStls.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <vsstls_pandas_table>`.
    
.. autoclass:: qupled.classic.VSStls
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.classic.VSStls.Input 
    :members:
    :inherited-members:
    :exclude-members: print, isEqual

.. _vsstls_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - alpha
     - ndarray
     - The free parameter
   * - fxcGrid
     - ndarray
     - The coupling parameter grid
   * - fxci
     - ndarray (2D)
     - The free energy integrand
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function
       
Hybrid schemes
--------------

ESA scheme
~~~~~~~~~~

The :obj:`qupled.classic.ESA` class is used to setup and perform all the necessary calculations
for the solution of the `Effective Static Approximation <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.ESA.Input`.
After the solution is completed the results are written to an hdf file in the form of :ref:`pandas dataframes <rpa_pandas_table>`.

.. autoclass:: qupled.classic.ESA
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.classic.ESA.Input
    :members:
    :show-inheritance:
		      
Quantum schemes
---------------

Qstls scheme
~~~~~~~~~~~~

The :obj:`qupled.quantum.Qstls` class is used to setup and perform all the necessary calculations
for the solution of the `Qstls scheme <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.quantum.Qstls.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <qstls_pandas_table>`.
     
.. autoclass:: qupled.quantum.Qstls
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.quantum.Qstls.Input 
    :members:
    :inherited-members:
    :exclude-members: print, isEqual

.. _qstls_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - adr
     - ndarray (2D)
     - Auxiliary density response
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function
       
Qstls-IET scheme
~~~~~~~~~~~~~~~~

The :obj:`qupled.quantum.QstlsIet` class is used to setup and perform all the necessary calculations
for the solution of the `Qstls-IET schemes <https://pubs.aip.org/aip/jcp/article/158/14/141102/
2877795/Quantum-version-of-the-integral-equation-theory>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.quantum.QstlsIet.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <qstlsiet_pandas_table>`.

.. autoclass:: qupled.quantum.QstlsIet
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.quantum.QstlsIet.Input 
    :members:
    :inherited-members:
    :exclude-members: print, isEqual

.. _qstlsiet_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - adr
     - ndarray (2D)
     - Auxiliary density response
   * - bf
     - ndarray
     - The bridge function adder
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function

QVSStls scheme
~~~~~~~~~~~~~~~~

The :obj:`qupled.quantum.QVSStls` class is used to setup and perform all the necessary calculations
for the solution of the QVSStls schemes.
The solution parameters are specified with a dedicated class called :obj:`qupled.quantum.QstlsIet.Input`.
After the solution is completed the results are written to an hdf file in the form of
:ref:`pandas dataframes <qvsstls_pandas_table>`.

.. autoclass:: qupled.quantum.QVSStls
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.quantum.QVSStls.Input 
    :members:
    :inherited-members:
    :exclude-members: print, isEqual

.. _qvsstls_pandas_table:
.. list-table:: Content of the pandas dataframe stored in the output file
   :widths: 25 25 50
   :header-rows: 1

   * - Item
     - Data Type
     - Description
   * - :ref:`info <stls_info_table>`
     - Pandas DataFrame
     - Information on the solution
   * - adr
     - ndarray (2D)
     - Auxiliary density response
   * - alpha
     - ndarray
     - The free parameter
   * - fxcGrid
     - ndarray
     - The coupling parameter grid
   * - fxci
     - ndarray (2D)
     - The free energy integrand
   * - idr
     - ndarray (2D)
     - Ideal density response
   * - sdr
     - ndarray
     - Static density response
   * - slfc
     - ndarray
     - Static local field correction
   * - ssf
     - ndarray
     - Static structure factor
   * - ssfHF
     - ndarray
     - Hartree-Fock static structure factor
   * - wvg
     - ndarray
     - Wave-vector grid
   * - rdf*
     - ndarray
     - Radial distribution function
   * - rdfGrid*
     - ndarray
     - Grid used to compute radial distribution function
       
Initial guess
-------------

.. autoclass:: qupled.classic.Stls.Guess
    :members:
    
.. autoclass:: qupled.quantum.Qstls.Guess
    :members:


Output
------

.. autoclass:: qupled.util.HDF
   :members:

