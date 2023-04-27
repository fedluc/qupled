Examples
========

The following examples illustrate some of the most common ways to use qupled
in order to solve a given dielectric scheme.

A simple STLS solution
----------------------

This example sets up a simple STLS calculation plots some of the results 
that are produced once the calculation are completed. In order not to clutter too
much the plot of the ideal density response we plot only the results for a few
matsubara frequencies. At the end of the example the data for the static structure
factor is extracted from the object used in the solution and printed on screen.

.. literalinclude:: ../examples/solveStls.py
   :language: python

Solving the  classical IET schemes
----------------------------------

This example shows how to solve two classical STLS-IET schemes. First the STLS-HNC
scheme is solved, then the properties of the solution object are changed and the
STLS-LCT scheme is sovled


.. literalinclude:: ../examples/solveStlsIet.py
   :language: python

Solving the quantum schemes
---------------------------

This example shows how to solve the quantum dielectric schemes QSTLS and QSTL-LCT.
Since these schemes can have a relatively high computational cost, in this example
we limit the number of matsubara frequencies to 16, we use 16 OMP threads to
speed up the calculation and we employ a segregated approach to solve the two-dimensional
integrals that appear in the schemes


.. literalinclude:: ../examples/solveQuantumSchemes.py
   :language: python

	      	      
