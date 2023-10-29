Examples
========

The following examples present some common use cases that show how to run qupled and how to post-process the results.

A simple STLS solution
----------------------

This example sets up a simple STLS calculation plots some of the results 
that are produced once the calculation are completed. In order not to clutter too 
much the plot of the ideal density response we plot only the results for a few 
matsubara frequencies. There are two ways to access the results of the calculation:
Directly from the object used to perform the calculation and from the output file
created at the end of the run. The example illustrates how the static structure factor
can be accessed with both these methods, other quantities can be accessed in the same
way.

.. literalinclude:: ../examples/solveStls.py
   :language: python

Solving the  classical IET schemes
----------------------------------

This example shows how to solve two classical STLS-IET schemes. First the STLS-HNC 
scheme is solved, then the properties of the solution object are changed and the 
STLS-LCT scheme is sovled

.. literalinclude:: ../examples/solveStlsIet.py
   :language: python

Solving the classical VS-STLS scheme
------------------------------------

This example shows how to solve the VS-STLS scheme for low coupling and zero
temperature. First the scheme is solved up to rs = 0.2, then the results are
plotted and then the calculation is resumed up to rs = 0.4. In the second
part of the calculation, the pre-computed value of the free energy integrand
available from the VS-STLS solution at rs = 0.2 is used in order to speed
up the calculation

.. literalinclude:: ../examples/solveVSStls.py
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

	      
Speed-up the quantum schemes
----------------------------

The calculations for the quantum schemes can be made significantly faster if part of the 
calculation of the auxiliary density response can be skipped. This can usually be done 
by passing in input the so-called 'fixed' component of the auxiliary density response. The 
fixed component of the auxiliary density response depends only on the degeneracy parameter 
and is printed to specific ouput files when a quantum scheme is solved. These output files can 
be used is successive calculations to avoid recomputing the fixed component and to speed-up 
the solution of the quantum schemes. The following two examples illustrate how this can be 
done for both the QSTLS and the QSTLS-IET schemes.

For the QSTLS scheme it is sufficient to pass a binary file containing the fixed component. 
This allows to obtain identical results (compare the internal energies printed at the end of 
the example) in a fraction of the time. We can also recycle the same fixed component for 
different coupling parameters provided that the degeneracy parameter stays the same. On the 
other hand, when changing the degeneracy parameter the fixed component of must also be upated 
otherwise the calculation fails as shown at the end of the example. 

.. literalinclude:: ../examples/fixedAdrQstls.py
   :language: python

For the QSTLS-IET schemes we must pass the name of two files: the binary file with the 
fixed auxiliary density response from the QSTLS scheme and a zip file containing a collection 
of binary files representing the fixed component for the QSTLS-IET scheme. Here the fixed 
component depends only on the degeneracy parameter but not on the coupling 
parameter and not on the theory used for the bridge function.

.. literalinclude:: ../examples/fixedAdrQstlsIet.py
   :language: python
