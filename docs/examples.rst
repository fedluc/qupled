Examples
========

The following examples present some common use cases that show how to run qupled and how to post-process the results.

Setup a scheme and analyze the output
-------------------------------------

This example sets up all the necessary objects to solve the RPA and ESA schemes and
shows how to access the information stored in the output files produced at the
end of the calculations

.. literalinclude:: ../examples/docs/solveRpaAndESA.py
   :language: python

A simple STLS solution
----------------------

This example sets up a simple STLS calculation and  plots some of the results 
that are produced once the calculation are completed. There are two ways to access
the results of the calculation: Directly from the object used to perform the calculation
or from the output file created at the end of the run. The example illustrates how
the static structure factor can be accessed with both these methods. Other quantities
can be accessed in the same way.

.. literalinclude:: ../examples/docs/solveStls.py
   :language: python

Solving the  classical IET schemes
----------------------------------

This example shows how to solve two classical STLS-IET schemes: the STLS-HNC and
the STLS-LCT schemes. The schemes are solved one after the other by simply
updating the properties of the solution object.

.. literalinclude:: ../examples/docs/solveStlsIet.py
   :language: python

Solving the classical VS-STLS scheme
------------------------------------

This example shows how to solve the VS-STLS scheme at finite temperature.
First the scheme is solved up to rs = 5.0, then the results are
plotted and then the calculation is resumed up to rs = 10. In the second
part of the calculation, the pre-computed value of the free energy integrand
available from the VS-STLS solution at rs = 5.0 is used in order to speed
up the calculation.

.. literalinclude:: ../examples/docs/solveVSStls.py
   :language: python

.. _solvingQuantumSchemes:

Solving the quantum schemes
---------------------------

This example shows how to solve the quantum dielectric schemes QSTLS and QSTLS-LCT. 
Since these schemes can have a relatively high computational cost, in this example 
we limit the number of matsubara frequencies to 16, we use 16 OMP threads to 
speed up the calculation and we employ a segregated approach to solve the two-dimensional 
integrals that appear in the schemes. 

.. literalinclude:: ../examples/docs/solveQuantumSchemes.py
   :language: python
	      	      
Define an initial guess
-----------------------

The following three examples show how to define an initial guess for the classical 
schemes (STLS and STLS-IET) and for the quantum schemes (QSTLS and QSTLS-IET). If 
an initial guess is not specified the code will use the default, namely zero static 
local field correction for the classical schemes and STLS static structure factor 
for the quantum schemes.

.. literalinclude:: ../examples/docs/initialGuessStls.py
   :language: python

In the following example we solve the QSTLS scheme twice and the second time we
specify the initial guess as the solution obtained from the first solution. Having
provided a nearly exact initial guess the scheme converges in a single iteration.
	      
.. literalinclude:: ../examples/docs/initialGuessQstls.py
   :language: python

The QSTLS-IET scheme  requires to specify an initial guess for the auxiliary density
response and the number of matsubara frequencies corresponding to such initial guess.
These specifications can be skipped in all other schemes.

.. literalinclude:: ../examples/docs/initialGuessQstlsIet.py
   :language: python

	      
Speed-up the quantum schemes
----------------------------

The quantum schemes can have a significant computational cost. There are two strategies
that can be employed to speed up the calculations:

* *Parallelization*: qupled supports both multithreaded calculations with OpenMP and
  multiprocessors computations with MPI. The number of OpenMP threads
  can be specified in input (as shown in :ref:`this example<solvingQuantumSchemes>`).
  Multiprocessor computations can be performed by running qupled as an MPI application:
  ``mpirun -n <number_of_cores> python3 <script_using_qupled>``. OpenMP and MPI can be
  used concurrently by setting both the number of threads and the number of cores.
 
* *Pre-computation*: The calculations for the quantum schemes can be made significantly
  faster if part of the calculation of the auxiliary density response can be skipped.
  This can usually be done by passing in input the so-called 'fixed' component of the
  auxiliary density response. The fixed component of the auxiliary density response depends
  only on the degeneracy parameter and is printed to specific ouput files when a quantum
  scheme is solved. These output files can  be used is successive calculations to avoid
  recomputing the fixed component and to speed-up  the solution of the quantum schemes.
  The following two examples illustrate how this can be done for both the QSTLS and
  the QSTLS-IET schemes.

.. literalinclude:: ../examples/docs/fixedAdrQstls.py
   :language: python

For the QSTLS scheme it is sufficient to pass a binary file containing the fixed component. 
This allows to obtain identical results (compare the internal energies printed at the end of 
the example) in a fraction of the time. We can also recycle the same fixed component for 
different coupling parameters provided that the degeneracy parameter stays the same. On the 
other hand, when changing the degeneracy parameter the fixed component of must also be upated 
otherwise the calculation fails as shown at the end of the example.

.. literalinclude:: ../examples/docs/fixedAdrQstlsIet.py
   :language: python

For the QSTLS-IET schemes we must pass the name of two files: the binary file with the 
fixed auxiliary density response from the QSTLS scheme and a zip file containing a collection 
of binary files representing the fixed component for the QSTLS-IET scheme. Here the fixed 
component depends only on the degeneracy parameter but not on the coupling 
parameter and not on the theory used for the bridge function.
