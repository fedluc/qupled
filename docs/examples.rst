Examples
========

The following examples present some common use cases that show how to run qupled and how to post-process the results.

Setup a scheme and analyze the output
-------------------------------------

This example solves the STLS and shows two ways in which to access the results: 
Directly from the object used to perform the calculation or from the database used
to store the results. 

.. literalinclude:: ../examples/docs/solve_stls.py
   :language: python

Define an initial guess
-----------------------

Qupled allows to specify a custom initial guess for any scheme.  This example shows how to define an 
initial guess for the STLS scheme, but the same approach can be used for any other scheme.

.. literalinclude:: ../examples/docs/initial_guess_stls.py
   :language: python

Read the results from the database
----------------------------------

Qupled stores the results of the calculations in a database that can be easily accessed to 
retrieve any quantity of interest. This example shows how to read the results from the database 
in order to plot the static structure factor for two different schemes.

.. literalinclude:: ../examples/docs/solve_rpa_and_esa.py
   :language: python

Post-process the results
------------------------

Apart from the implmented schemes, qupled also includes a set of post-processing tools that can be 
used to further analyze the results of the calculations and extract valuable information about the system.

Compute additional quantities from an existing object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how the post-processing tools can be used to compute the radial distribution function 
and the imaginary-time correlation function, two quantities that are not directly computed in the 
iterative calculations but that can be easily obtained from the results of the calculation.

.. literalinclude:: ../examples/docs/post_processing.py
   :language: python

Compute additional quantities from the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The post-processing tools can also be used to compute additional quantities by reading the necessary data 
from the database. This is particularly useful when revisiting previous calculations that were stored without
computing the additional quantities. This example shows how to compute the  radial distribution function
starting from run results stored in the database.

.. literalinclude:: ../examples/docs/correlation_functions.py
   :language: python

Finite Size Correction
~~~~~~~~~~~~~~~~~~~~~~

Qupled also includes a module that can be used to compute the finite size correction to the free energy. 
This example shows how to use this module to compute the finite size correction 
for a 3D system at :math:`r_s=5` and :math:`\Theta=1`.

.. literalinclude:: ../examples/docs/finite_size_correction.py
   :language: python

Computationally intestive schemes
---------------------------------

The quantum and the VS schemes are computationally more demanding than the classical schemes. 
The following examples show how what can be done to speed up the calculations for these schemes.

VS Schemes: Pre-compute the free energy integrand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates how to solve the classical VS-STLS scheme at finite temperature.
The calculation is first carried out up to :math:`r_s=2` and then resumed up to :math:`r_s=5`
reusing the pre-computed free energy integrand. VS-type schemes can be numerically demanding, so 
it is often convenient to be able to restart from a known state point using 
previously computed quantities.

.. literalinclude:: ../examples/docs/solve_vsstls.py
   :language: python

Quantum schemes: parallelization and pre-computation 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two strategies that can be employed to speed up the calculations of quantum schemes:

* *Parallelization*: qupled supports both multithreaded calculations with OpenMP and
  multiprocessors computations with MPI. OpenMP and MPI can be
  used concurrently by setting both the number of `threads` and the number of `processes` in the 
  input dataclasses. The following example shows how to solve the quantum schemes using OpenMP parallelization.

.. literalinclude:: ../examples/docs/solve_quantum_schemes.py
   :language: python 
 
* *Pre-computation*: The calculations for the quantum schemes can be made significantly
  faster if part of the calculation of the auxiliary density response can be skipped.
  Qupled will look into the database used to store the results to try to find the 
  necessary data to skip the full calculation of the auxiliary density response. Try 
  to run the following example and notice how the second calculation is much faster
  than the first one.

.. literalinclude:: ../examples/docs/fixed_adr.py
   :language: python
