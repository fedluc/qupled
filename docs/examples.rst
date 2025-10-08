Examples
========

The following examples present some common use cases that show how to run qupled and how to post-process the results.

Setup a scheme and analyze the output
-------------------------------------

This example sets up all the necessary objects to solve the RPA and ESA schemes and
shows how to access the information stored in the output files produced at the
end of the calculations

.. literalinclude:: ../examples/docs/solve_rpa_and_esa.py
   :language: python

Different ways to access the results
------------------------------------

This example solves the STLS and shows two ways in which to access the results: 
Directly from the object used to perform the calculation or from the database used
to store the results. 

.. literalinclude:: ../examples/docs/solve_stls.py
   :language: python

Solving the VS schemes
----------------------

This example demonstrates how to solve the classical VS-STLS scheme at finite temperature.
The calculation is first carried out up to :math:`r_s=2` and then resumed up to :math:`r_s=5`
reusing the pre-computed free energy integrand.
Because VS-type schemes can be numerically demanding, it is often convenient to be able to 
restart from a known state point using previously computed quantities. This approach not only 
saves computational time but also improves the robustness of the solution process.

.. literalinclude:: ../examples/docs/solve_vsstls.py
   :language: python
         
Define an initial guess
-----------------------

Qupled allows to specify a custom initial guess for any scheme.  This exampl shows how to define an 
initial guess for the STLS scheme, but the same approach can be used for any other scheme.

.. literalinclude:: ../examples/docs/initial_guess_stls.py
   :language: python
         
Speed-up the quantum schemes
----------------------------

The quantum schemes can have a significant computational cost. There are two strategies
that can be employed to speed up the calculations:

* *Parallelization*: qupled supports both multithreaded calculations with OpenMP and
  multiprocessors computations with MPI. OpenMP and MPI can be
  used concurrently by setting both the number of threads and the number of cores in the 
  input dataclasses. Use `threads` to set the number of OMP threads and `processes` to
  set the number of MPI processes, as shown in the following example

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