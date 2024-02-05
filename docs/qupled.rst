The package
===========

The package is written in both Python and C++. Python can be
used to setup, run, save and postprocess a dielectric scheme of choice. C++
is used to perform the actual calculations involved in the solution of the schemes.
Communication between Python and C++ is handled via
`python::boost <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_.
The following sections illustrate the classes used to solve classical and
quantum schemes and the classes used to handle the inputs, initial guesses and outputs.

Classic schemes
---------------

.. autoclass:: qupled.classic.Rpa
    :members:
    :undoc-members:

.. autoclass:: qupled.classic.Stls
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: qupled.classic.StlsIet
    :members:
    :undoc-members:
    :show-inheritance:
       
.. autoclass:: qupled.classic.VSStls
    :members:
    :show-inheritance:

Hybrid schemes
--------------

.. autoclass:: qupled.classic.ESA
    :members:
    :show-inheritance:
    
Quantum schemes
---------------
       
.. autoclass:: qupled.quantum.Qstls
    :members:
    :show-inheritance:

.. autoclass:: qupled.quantum.QstlsIet
    :members:
    :show-inheritance:
				
Input
-----

.. autoclass:: qupled.qupled.RpaInput
    :members:
    :undoc-members:

.. autoclass:: qupled.qupled.StlsInput
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: qupled.qupled.VSStlsInput
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: qupled.qupled.FreeEnergyIntegrand
    :members:
    :undoc-members:

.. autoclass:: qupled.qupled.QstlsInput
    :members:
    :undoc-members:
    :show-inheritance:

Initial guess
-------------

.. autoclass:: qupled.qupled.SlfcGuess
    :members:
    :undoc-members:
    
.. autoclass:: qupled.qupled.QstlsGuess
    :members:
    :undoc-members:


Output
------

.. autoclass:: qupled.util.Hdf
   :members:
