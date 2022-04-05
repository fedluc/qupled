# Dielectric formalism for strongly coupled plasmas

## Introduction

STLS can be used to compute the properties of quantum one component plasmas via theoretical approaches based on the dielectric formalism. The theoretical approaches which can be solved with STLS include:

* The classical STLS scheme as discussed by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278)
* The classical VS-STLS scheme discussed by [Vashishta and Singwi](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875) extended to finite temperatures
* The classical STLS-HNC scheme as discussed by [Tanaka](https://aip.scitation.org/doi/full/10.1063/1.4969071)
* The classical STLS-IET scheme as discussed by [Tolias and collaborators](https://aip.scitation.org/doi/10.1063/5.0065988)
* The quantum STLS (qSTLS) scheme as discussed by [Schweng and Böhm](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037)
* The quantum qSTLS-IET scheme 
 
## Limitations

Ground state (zero temperature) calculations are available only for the classical schemes (STLS, STLS-HNC, STLS-IET and VS-STLS).

## Units

All the calculations are performed in normalized units. The wave vectors are normalized to the Fermi wave-vector and the frequencies are normalized to 2&pi;E<sub>F</sub>/h. Here E<sub>F</sub> is the Fermi energy and h is Planck's constant.

## Compiling

The code can be compiled with gcc and with the [make file](Makefile) provided with the source code. In order to correctly compile the code it is necessary that the following libraries are installed

* [GNU scientific library](https://www.gnu.org/software/gsl/). This library must be explicitly included by calling `make GSL="PATH"`, where `PATH` is the path to the folder containing the header files of the library. Alternatively, it is possible to change the default value of the GSL variable in the [make file](Makefile) and then to simply compile via `make`.
* [OpenMP library](https://en.wikipedia.org/wiki/OpenMP). In most cases it is not necessary to explicitly include the path to this library. However, if this should be necessary, it should be possible to do so by modifying the `INCLUDE` variables which appear in the [make file](Makefile).

## Running

SLTS can run in three different modes: `static`,`dynamic` and `restart`. Each of these modes is described in detail below and `static` is the default working mode of the code. The user is free to specify how the code should run via the run-time option `mode`, but only one mode at the time can be executed. All the integrals that appear in the various dielectric schemes are computed with the doubly-adaptive Clenshaw-Curtis quadrature scheme as implemented in the [CQUAD](https://www.gnu.org/software/gsl/doc/html/integration.html) function of the GSL library. The run-time options that can be specified by the user are summarized in a [dedicated section](#run-time-options).

### Static mode 

In the static mode, SLTS computes the static structure factor and related static and thermodynamic properties with an iterative procedure that depends on the theory specified by the user:

#### Classical schemes (STLS, STLS-HNC, STLS-IET)
  * The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method. The initial guess for the bisection method is controlled via the input option `mu-guess`
  * The ideal density response is computed for various values of the wave-vector and of the Matsubara frequency
  * The static structure factor and the static local field correction are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||G<sub>i</sub>(x) - G<sub>i-1</sub>(x)|| < &epsilon; is satisfied between two successive iterations. Here G(x) is the static local field correction and &epsilon; is a tolerance specified in input via the option `min-err`. Unless specified otherwise (see option `stls-guess`), the initial guess for the iterative calculations corresponds to the RPA solution, G(x) = 0

#### Classical VS-STLS scheme
  * The free parameter used to enforce the compressibility sum-rule (&alpha;) is read from input (see option `vs-alpha`)
  * The static structure factor and the static local field correction are computed with the an iterative procedure similar to the one employed for the other classical schemes (see above). In order to compute the state point derivatives that appear in the expression for the static local field correction, nine state points are solved simultaneously
  * The interaction energy is computed from the static structure factor
  * The exchange free energy is computed from the internal energy
  * A new value for the free parameter used to enforce the compressibility sum-rule is obtained from the interaction energy and from the exchange free energy. If the condition |&alpha;<sub>i</sub>(x) - &alpha;<sub>i-1</sub>(x)|/ &alpha;<sub>i</sub>(x) < &epsilon; is satisfied, the calculation is complete. Otherwise, the new value for &alpha; is used to start a new computational cycle. The tolerance &epsilon; used to check the convergence of the free parameter is specified via the option `vs-min-err` and does not have to be equal to the tolerance used to compute the structural properties

#### Quantum schemes (qSTLS, qSTLS-IET)
  * The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method
  * The ideal density response is computed for various values of the wave-vector and of the Matsubara frequency 
  * The fixed component of the auxiliary density response (defined as the product between the dynamic local field correction and the ideal density response) that does not depend explicitly on the static structure factor is computed and stored. 
  * The static structure factor and the part of the auxiliary density response that depends explicitly on the static structure factor are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||S<sub>i</sub>(x) - S<sub>i-1</sub>(x)|| < &epsilon; is satisfied between two successive iterations. Here S(x) is the static structure factor and &epsilon; is a tolerance specified via the option `min-err`. Unless specified otherwise (see option `qstls-guess`), the initial guess for the iterative calculations corresponds to zero auxiliary density response.

#### Output

Once the iterative procedure is completed, the results are written to a set of dedicated output files which include (n<sub>x</sub> is the number of wave-vectors used in the solution and n<sub>l</sub> is the number of Matsubara frequencies):
  
  * ssf_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and two columns. The first column is the wave vector, the second column is the static structure factor
  * ssfHF_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and two columns. The first column is the wave vector, the second column is the static structure factor within the Hartree-Fock approximation
  * slfc_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and two columns. The first column is the wave vector, the second column is the static local field correction
  * sdr_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and two columns. The first column is the wave vector, the second column is the static density response
  * idr_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and n<sub>l</sub> columns. Each column corresponds to the ideal density response for a given value of the Matsubara frequency, starting from l = 0 in the first column
  * rdf_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and two columns. The first column is the inter-particle distance (with a fixed resolution of 0.01), the second column is the radial distribution function
  * uint_rs\*\_theta\*\*\_\*\*\*.dat, a text file with three elements: the quantum coupling parameter, the quantum degeneracy parameter and the interaction energy.
  * restart_rs\*\_theta\*\*\_\*\*\*.bin, a binary file which can be used to restart simulations that were interrupted before they reached convergence or that can be employed to provide an initial guess for subsequent calculations. In order to understand how to use this binary file, see the options `--stls-restart` and `--qstls-restart`
   
 For the VS-STLS, additional files are produced: 
 
 * alpha_csr_rs\*\_theta\*\*\_\*\*\*.dat, a text file with three elements: the quantum coupling parameter, the quantum degeneracy parameter and the free parameter used to enforce the compressibility sum rule.
  * fxc_rs\*\_theta\*\*\_\*\*\*.dat, a text file with three elements: the quantum coupling parameter, the quantum degeneracy parameter and the exchange free energy.
  * thermo_int_rs\*\_theta\*\*\_\*\*\*.bin, a binary file with the free energy integrand that can be used for subsequent calculations via the option `vs-thermo-file` 
  
 For the quantum schemes, additional files are produced: 
 
  * adr_rs\*\_theta\*\*\_\*\*\*.dat, a text file with n<sub>x</sub> rows and n<sub>l</sub> columns. Each column corresponds to the auxiliary density response for a given value of the Matsubara frequency, starting from l = 0 in the first column
  * fixed_rs\*\_theta\*\*\_\*\*\*.bin, a binary file which can be used to avoid recomputing the fixed component of the auxiliary density response in the qSTLS scheme via the option `--qstls-fix`  () **(only for the qSLTLS scheme)** 
  * adr_iet_fixed_theta\*\*\_xx\*\*\*\*.bin, n<sub>x</sub> binary files which can be used to avoid recomputing the fixed component of the auxiliary density response in the qSTLS-IET scheme via the option `--qstls-iet-fix` **(only for the qSTLS-IET scheme)**
  
In the above, \* corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter, \*\*\* to the dielectric scheme that was solved and \*\*\*\* is the wave-vector value.

### Dynamic mode

In the dynamic mode, SLTS computes the dynamic structure factor and related properties for a given dielectric scheme. The necessary structural input is passed via an input file specified via the option `--dyn-struct`. For the classical schemes, the structural input consists of the static local field correction. For the quantum schemes, the structural input consists of the static structure factor. The format of the files must be consistent with the one described [here](#output).

The dynamic properties are computed over a grid of frequencies whose cutoffs and resolution are controlled via the run-time options `dyn-Wmin`, `dyn-Wmax` and `dyn-dW`, and for all the wave-vectors contained in the input file with the structural properties. 

#### Output

Once the calculation is completed, the results are written to a set of dedicated output files which include (n<sub>W</sub> is the number of frequencies used in the solution):
 
  * dsf_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a text file with n<sub>W</sub> rows and two columns. The first column is the frequency, the second column is the dynamic structure factor for the wave-vector specified via the run-time option `dyn-xtarget`
  * isf_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a text file with n<sub>W</sub> rows and two columns. The first column is the imaginary time (defined between 0 and 1), the second column is the intermediate scattering function for the wave-vector specified via the run-time option `dyn-xtarget`
  * idr_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a text file with n<sub>W</sub> rows and three columns. The first column is the frequency, the second column the real part of the ideal density response and the third column the imaginary part of the ideal density response. All the density responses refer to the wave-vector specified via the run-time option `dyn-xtarget`
  * dynamic_restart_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a binary file which contains the ideal density response for all the wave-vectors contained in the input file with the structural properties. This file can be used in subsequent calculations with different values of `dyn-xtarget` in order to compute the dynamic properties at different wave vectors without having to recompute the density response from scratch. In order to understand how to use this binary file, see the option `--dyn-restart`

 For the quantum schemes, additional files are produced: 
 
  * adr_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a text file with n<sub>W</sub> rows and three columns. The first column is the frequency, the second column the real part of the auxiliary density response and the third column the imaginary part of the auxiliary density response. All the density responses refer to the wave-vector specified via the run-time option `dyn-xtarget`. 
  
  * dynamic_restart_rs\*\_theta\*\*\_x\*\*\*\_\*\*\*\*\.dat, a binary file which contains the ideal and auxiliary density response for all the wave-vectors contained in the input file with the structural properties. This file can be used in subsequent calculations with different values of `dyn-xtarget` in order to compute the dynamic properties at different wave vectors without having to recompute the density response from scratch. In order to understand how to use this binary file, see the option `--dyn-restart`

In the above, \* corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter,  \*\*\* to the wave-vector value specified via  `dyn-xtarget` and \*\*\*\* to the dielectric scheme that was solved.

### Restart mode

Upon completion of a calculation in the static mode, STLS writes a binary restart file that can be used to restart the simulation (or as a guess file for subsquent calculations). However, should this restart file not be available, the restart mode can be used to create the missing binary file. In the restart mode the information necessary to construct the restart file is read from two text files provided by the user via the option `--restart-files`.  The format of the text files depends on the type of restart file that should be generated: 

 * Restart file for the classical schemes: The first text file must contain the static structure factor, the second file must contain the static local field correction. The format of the files must be consistent with the one described [here](#output)

* Restart file for the quantum schemes: The first text file must contain the static structure factor, the second file must contain the auxiliary density response. The format of the files must be consistent with the one described [here](#output)

The information concerning the grid and number of Matsubara frequencies is inferred from the structure of the text files, but information concerning the state point and the scheme is obtained from the input parameters specified by the user. A summary of the parameters used to construct the binary file is written on screen. After reading the text files, STLS writes a binary file called restart_rs\*\_theta\*\*\_\*\*\*.bin, where \*  corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter and \*\*\* to the dielectric scheme  specified in input. 

## Run-time options

The following command line options can be employed to control the calculations performed by STLS (the same information can also be retrieved by running STLS with the option `--help`) :

  * `--debug-input` can be used to print the content of the entire input structure to the screen, mainly useful for debugging purposes. The entire input structure is printed if `--debug-input=1` is specified, otherwise only a summary of the most relevant input parameters is printed. Default `--debug-input=0`
  
  * `--dx` specifies the  resolution for wave-vector grid. Default `--dx=0.1`
   
  * `--dyn-dw` specifies the resolution for the frequency grid for the dynamic properties. Default `--dyn-dw=0.1`'

  * `--dyn-restart` specifies the name of the binary file used to load the density response (ideal and auxiliary) for the calculation of the dynamic properties. Default: no file is specified and the density response is computed from scratch

  * `--dyn-struct` specifies the name of the text file used to load the structural properties for the calculation of the dynamic properties. Default: for the classical schemes the static local field correction file slfc_rs\*\_theta\*\*\_\*\*\*.dat is loaded, for the quantum schemes the static structure factor file ssf_rs\*\_theta\*\*\_\*\*\*.dat is loaded. Here,  \* corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter and \*\*\* to the dielectric scheme.

  * `--dyn-wmax` specifies the upper cutoff for the frequency grid used in the calculation of the dynamic properties. Default `--dyn-wmax=20.0`

  * `--dyn-wmin` specifies the lower cutoff for the frequency grid used in the calculation of the dynamic properties. Default `--dyn-wmin=0.0`
  
  * `--dyn-xtarget` specifies the wave-vector used to output the dynamic properties, only one value can be specified. Default `--dyn-xtarget=1.0`
  
  * `--iet-mapping` specifies what mapping to use between classical state points specified by &Gamma; and quantum state points specified by (r<sub>s</sub>, &theta;) for the IET-based theories (STLS-IET, qSTLS-IET). Default `--iet-mapping=standard`. Accepted options are: 
     * `standard` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/&theta;  (NOTE: this cannot be used in the ground state)
     * `sqrt` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/(1 + &theta;<sup>2</sup>)<sup>1/2</sup>
     * `linear` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/(1 + &theta;)

  * `--iter`  specifies the maximum number of iterations to employ in the iterative procedure used to compute the static structure factor.  Default `--iter=1000`.
  
  * `--min-err` specifies the minimum error for convergence in the iterative procedure used to compute the static structure factor.  Default `--min-err=1e-5`.
  
  * `--mix` specifies the mixing parameter for iterative solution.  Default `--mix=0.1`

  * `--mode` specifies the working mode of the code. Default: `--mode=static`

  * `--mu-guess` specifies the initial guess for the bisection procedure used to compute the chemical potential.  Default `--mgu-guess=-10,10`
  
  * `--nl`  specifies the number of Matsubara frequencies.  Default `--nl=128`

  * `--omp` specifies how many omp threads to use for the parts of the code that can be executed in parallel. Since the classical schemes can usually be run efficiently with a single thread, parallelization is available only for the quantum schemes. Default `omp=1`

  * `--qstls-fix` specifies the name of the binary file used to load the fixed component of the auxiliary density response function in the solution of the qSTLS scheme. If a file name is not specified, the fixed component of the auxiliary density response is computed from scratch and a significant increase in the computational cost for the solution of the qSTLS scheme can be expected. It is useful to note that state points with the same degeneracy parameter share the same fixed component of the auxiliary density response. Default: no file name is specified

* `--qstls-iet-fix` specifies the name of the folder containing the binary files used to load the fixed component of the auxiliary density response function in the solution of the qSTLS-IET scheme.  The files are assumed to have the name psi_fixed_theta\*_xx\**_.bin, where \* is the value of the quantum degeneracy parameter and \** is the value of the wave-vector. One file for each wave-vector is expected. If a folder name is not specified, the fixed component of the auxiliary density response is computed from scratch and a significant increase in the computational cost for the solution of the qSTLS-IET scheme can be expected. It is useful to note that state points with the same degeneracy parameter share the same fixed component of the auxiliary density response. Default: no file name is specified

* `--qstls-iet-static` specifies how the auxiliary density response should be computed within the qSTLS-IET scheme. If `qstls-iet-static=0` the auxiliary density response is computed within the fully dynamic approximation. Otherwise, If `qstls-iet-static=1` the auxiliary density response is computed within the partially dynamic approximation. Default `qstls-iet-static=0`

  * `--qstls-restart` specifies the name of the binary file used as initial guess for the solution of the quantum schemes. If a file name is not specified, the iterative scheme is started by setting the auxiliary density response to zero (which is equivalent to set the dynamic local field correction to zero). Default: no file name is specified
 
 * `--restart-files` specifies the names of the two text files used to construct the binary files when the code is run in restart mode. Default: no file names are specified and the code is run in static mode
  
  * `--rs`  specifies the  quantum coupling parameter. Default `--rs=1.0`
  
  * `--stls-restart` specifies the name of the binary file used as initial guess for the solution of the classical schemes. If a file name is not specified, the iterative scheme is started with the RPA approximation (zero static local field correction). Default: no file is specified 
  
  * `--theory` specifies which dielectric scheme should be solved . Default `--theory=STLS`. Accepted options are: 
     * `STLS` (classical STLS scheme), 
     * `STLS-HNC` (classical STLS-HNC scheme) 
     * `STLS-IET-IOI` (classical STLS-IET scheme with the [IOI bridge function](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.1051))
     * `STLS-IET-LCT` (classical STLS-IET scheme with the [LCT bridge function](https://arxiv.org/abs/2108.09574))
     * `QSTLS` (quantum STLS scheme). 
     * `QSTLS-IET-HNC` (quantum STLS-HNC scheme). 
     * `QSTLS-IET-IOI` (quantum STLS-IET scheme with the [IOI bridge function](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.1051))
     * `QSTLS-IET-LCT` (quantum STLS-IET scheme with the [LCT bridge function](https://arxiv.org/abs/2108.09574))
     
  * `--Theta` specifies the  quantum degeneracy parameter. Default `--Theta=1.0`
  
  * `--vs-alpha` specifies the initial guess for the free parameter used to enforce the compressibility sum-rule in the VS-STLS scheme. Default `--vs-alpha=0.5`
  
  * `--vs-drs` speficies the resolution used to compute the derivatives with respect to the coupling parameter and the exchange free energy integrand in the VS-STLS scheme. Default `--vs-drs=0.01`
  
  * `--vs-dt` speficies the resolution used to compute the derivatives with respect to the degeneracy parameter in the VS-STLS scheme. Default `--vs-dt=0.01`
  
  * `--vs-min-err` speficies the minimum error for convergence in the iterations for the free parameter used to enforce the compressibility sum-rule in the VS-STLS scheme. Default `--vs-alpha=1e-3`
    
  * `--vs-mix` speficies the mixing parameter in the iterations for the free parameter used to enforce the compressibility sum-rule in the VS-STLS scheme . Default `--vs-mix=1.0`  
  
  * `--vs-solve-csr` speficies whether to enforce the compressibility sum-rule in the VS-STLS scheme or not. If this parameter is set to 0, the self consistent calculation for the free parameter in the VS-STLS is by-passed completely and the structural properties are determined via the free parameter specified with `vs-alpha`. Default `--vs-solve-csr=1`
  
  * `--vs-thermo` speficies the name of the binary file used to load part of the exchange free energy integrand in the VS-STLS scheme. The binary file with the free energy integrand is written at the end of any successfull VS-STLS calculation, see [here](#output). If no file name is given, the free energy integrand is computed starting from 0. Note that the same value for the free parameter used to enforce the compressibility sum-rule, &alpha;, is adopted for all the state points that are not included in the imported binary file. Hence, if &alpha; is expected to vary significantly with the coupling parameter, computing the free energy integrand from 0 could lead to erroneous results.  Default: no file name is specified.

  * `--xmax` specifies the cutoff for wave-vector grid. Default `--xcut=20`
 
  ## Constants
  
  The file [utils.h](https://github.com/fedluc/STLS/blob/master/utils.h) contains five constants defined at compile time which are used to define the relative errors and maximum number of iterations for the root solvers and quadrature schemes. For state points where it is hard to obtain convergence it can be useful to adjust the values of such constants.
 

