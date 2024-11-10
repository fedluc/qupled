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
After the solution is completed the results are written to an hdf file in the form of.

.. autoclass:: qupled.classic.Rpa
    :members:
    :inherited-members:
    :exclude-members: Input, rdf

.. autoclass:: qupled.qupled.Rpa
    :members:
