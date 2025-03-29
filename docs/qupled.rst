The package
===========

The qupled package is implemented using both Python and C++. Python is responsible for setting up, 
running, and storing the results of a chosen dielectric scheme in an SQL database. C++ handles 
the computational tasks required to solve the schemes. Communication between Python and C++ is 
facilitated through the `Boost <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_ 
library. The following sections provide an overview of the classes used for solving classical 
and quantum schemes, as well as those for managing inputs, initial guesses, and outputs.


Classic schemes
---------------

Rpa scheme
~~~~~~~~~~

The :obj:`qupled.classic.Rpa` class is used to setup and perform all the necessary calculations
for the solution of the `Random-Phase Approximation <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.Rpa.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.classic.Rpa.Result` 
and written to the output database.

.. autoclass:: qupled.classic.Rpa
    :members:
    :inherited-members:
    :exclude-members: Input, Result

.. autoclass:: qupled.classic.Rpa.Input
    :inherited-members:
    :members:
    :exclude-members: print, isEqual, to_native

.. autoclass:: qupled.classic.Rpa.Result
    :inherited-members:
    :members:
    :exclude-members: compute_rdf, from_native


Stls scheme
~~~~~~~~~~~		      

The :obj:`qupled.classic.Stls` class is used to setup and perform all the necessary calculations
for the solution of the `Stls scheme <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.Stls.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.classic.Stls.Results` 
and written to the output database.

.. autoclass:: qupled.classic.Stls
    :members:
    :inherited-members:
    :exclude-members: Input, Result

.. autoclass:: qupled.classic.Stls.Input	       
    :show-inheritance:
    :members:
    :exclude-members: print, isEqual, to_native
      
.. autoclass:: qupled.classic.Stls.Result
    :show-inheritance:
    :members:
    :exclude-members: compute_rdf, from_native


Stls-IET schemes
~~~~~~~~~~~~~~~~	      

The :obj:`qupled.classic.StlsIet` class is used to setup and perform all the necessary calculations
for the solution of the `Stls-IET schemes <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
The solution parameters are specified with a dedicated class called :obj:`qupled.classic.StlsIet.Input`.
After the solution is completed the results are stored in an object :obj:`qupled.classic.StlsIet.Results` 
and written to the output database.
    
.. autoclass:: qupled.classic.StlsIet
    :members:
    :inherited-members:
    :exclude-members: Input, Result

.. autoclass:: qupled.classic.StlsIet.Input 
    :show-inheritance:
    :members:
    :exclude-members: print, isEqual, to_native

.. autoclass:: qupled.classic.StlsIet.Result
    :show-inheritance:
    :members:
    :exclude-members: compute_rdf, from_native
