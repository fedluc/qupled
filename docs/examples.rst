Examples
========

The following examples present some common use cases that show how to run qupled and how to post-process the results.

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

	      	      
Define an initial guess
-----------------------

The following three examples show how to define an initial guess for the classical 
schemes (STLS and STLS-IET) and for the quantum schemes (QSTLS and QSTLS-IET). If 
an initial guess is not specified the code will use the default, namely zero static 
local field correction for the classical schemes and STLS static structure factor 
for the quantum schemes.

.. literalinclude:: ../examples/initialGuessStls.py
   :language: python

Here it's important to notice how the initial guess in the quantum schemes is 
assigned to the ``qInput`` field of the object used to describe the dielectric 
scheme. 
	   
.. literalinclude:: ../examples/initialGuessQstls.py
   :language: python

One should also pay attention to the fact the QSTLS-IET scheme 
requires to specify an initial guess for the auxiliary density response and the number
of matsubara frequencies corresponding to such initial guess. These specifications 
can be skipped in all other schemes.

.. literalinclude:: ../examples/initialGuessQstlsIet.py
   :language: python
