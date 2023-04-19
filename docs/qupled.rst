The qupled package
==================

The package is written in both python and C++. :ref:`Python classes` can be
used to setup, run, save and postprocess a dielectric scheme of choice. C++
is used to perform the actual calculations involved in the solution of the scheme.
The user should not worry about the details of the C++ code except for a few
classes that are exposed via
`python::boost <https://www.boost.org/doc/libs/1_80_0/libs/python/doc/html/index.html>`_
and that are discussed in the section :ref:`Exposed C++ methods` 

.. _Python classes:

Python classes 
--------------

.. autoclass:: qupled.Static.Stls
    :members:
    :undoc-members:

.. autoclass:: qupled.Static.StlsIet
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: qupled.Static.Qstls
    :members:
    :show-inheritance:

.. autoclass:: qupled.Static.QstlsIet
    :members:
    :show-inheritance:

.. _Exposed C++ methods:

Exposed C++ methods 
-------------------

.. autoclass:: qupled.qupled.Input
    :members:
    :undoc-members:

.. autoclass:: qupled.qupled.SlfcGuess
    :members:
    :undoc-members:

.. autoclass:: qupled.qupled.StlsInput
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: qupled.qupled.QstlsGuess
    :members:
    :undoc-members:
    
.. autoclass:: qupled.qupled.QstlsInput
    :members:
    :undoc-members:
    
.. autoclass:: qupled.qupled.Stls
    :members:
    :undoc-members:
       
.. autoclass:: qupled.qupled.Qstls
    :members:
    :undoc-members:
    :show-inheritance:

.. automethod:: qupled.qupled.computeRdf

.. automethod:: qupled.qupled.computeInternalEnergy
