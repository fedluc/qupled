# Dielectric formalism for strongly coupled plasmas

## Introduction

STLS can be used to compute the static and thermodynamic properties of quantum one component plasmas via theoretical approaches based on the dielectric formalism. The theoretical approaches which can be solved with STLS include:

* The classical STLS scheme as discussed by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278)
* The classical STLS-HNC scheme as discussed by [Tanaka](https://aip.scitation.org/doi/full/10.1063/1.4969071)
* The classical STLS-IET scheme as discussed by [Tolias and collaborators](https://aip.scitation.org/doi/10.1063/5.0065988)
* The quantum STLS (qSTLS) scheme as discussed by [Schweng and Böhm](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037)
* The quantum qSTLS-IET approach 
 
## Limitations

STLS can only be employed for finite-temperature systems. Calculations of ground state (zero temperature) properties will be implemented in the future.

## Compiling

The code can be compiled with gcc and with the [make file](Makefile) provided in the source directory. Please note that in order to correctly compile the code it is necessary that the following libraries are installed

* [GNU scientific library](https://www.gnu.org/software/gsl/). This library must be explicitly included by calling `make GSL="PATH"`, where `PATH` is the path to the folder containing the header files of the library. Alternatively, it is possible to change the default value of the GSL variable in the Makefile and then to simply compile via `make`.
* [OpenMP library](https://en.wikipedia.org/wiki/OpenMP). In most cases it is not necessary to explicitly include the path to this library. However, if this should be necessary, it should be possible to do so by modifying the `INCLUDE` variables which appear in the Makefile.

## Running 

Given a state point defined via the quantum degeneracy parameter (Theta) and via the quantum coupling parameter (r<sub>s</sub>), SLTS computes the static structure factor with the following procedure:

* For the classical schemes (STLS, STLS-HNC, STLS-IET)
  * The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method
  * The ideal density response is computed for various values of the wave-vector and of the Matsubara frequency
  * The static structure factor and the static local field correction are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||G<sub>i</sub>(x) - G<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here G(x) is the static local field correction and epsilon is a tolerance specified in input


* For the quantum schemes (qSTLS, qSTLS-IET)
  * The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method
  * The ideal density response is computed for various values of the wave-vector and of the Matsubara frequency 
  * The fixed component of the auxiliary density response (defined as the product between the dynamic local field correction and the ideal density response) that does not depend explicitly on the static structure factor is computed and stored. 
  * The static structure factor and the part of the auxiliary density response that depends explicitly on the static structure factor are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||S<sub>i</sub>(x) - S<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here S(x) is the static structure factor and epsilon is a tolerance specified in input.

All the integrals are computed with the doubly-adaptive Clenshaw-Curtis quadrature scheme as implemented in the [CQUAD](https://www.gnu.org/software/gsl/doc/html/integration.html) function of the GSL library. Once the iterative procedure is completed, the results are written to a set of dedicated output files (see the section "Output")

The following command line options can be employed to control the calculations performed by STLS (the same information can also be retrieved by running STLS with the option `--help`) :

  * `--debug-input` can be used to print the content of the entire input structure to the screen, mainly useful for debugging purposes. The entire input structure is printed if `--debug-input=1` is specified, otherwise only a summary of the most relevant input parameters is printed. Default `--debug-input=0`
  
  * `--dx` specifies the  resolution for wave-vector grid. Default `--dx=0.01`
  
  * `--guess-files` specifies the names of the two text files used to construct the binary files that can be supplied as an initial guess for the code. More information on the files for the initial guess are given in the section **Guess and restart**. Default `--guess-files=NO_FILE,NO_FILE` (no files are specified)
  
  * `--guess-write` can be used to run the code in "guess" mode by setting `--guess-write=1`. More information on what this means is given in the section **Guess and restart**. Default `--guess-write=0` (the code runs in the normal mode and solves the theory specified by `--theory`)
  
  * `--iet-mapping` specifies what mapping to use between classical state points specified by &Gamma; and quantum state points specified by (r<sub>s</sub>, &theta;) for the IET-based theories (STLS-IET, qSTLS-IET). Accepted options  are: 
     * `standard` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/&theta;
     * `sqrt` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/(1 + &theta;<sup>2</sup>)<sup>1/2</sup>
     * `linear` : &Gamma; = 2&lambda;<sup>2</sup>r<sub>s</sub>/(1 + &theta;)

  * `--iter`  specifies the maximum number of iterations to employ in the iterative procedure used to compute the static structure factor.  Default `--iter=1000`
  
  * `--min-err` specifies the minimum error for convergence in the iterative procedure used to compute the static structure factor.  Default `--min-err=1e-5`
  
  * `--mix` specifies the mixing parameter for iterative solution.  Default `--mix=0.1`

  * `--mu-guess` specifies the initial guess for the bisection procedure used to compute the chemical potential.  Default `--mgu-guess=-10,10`
  
  * `--nl`  specifies the number of Matsubara frequencies.  Default `--nl=128`.

  * `--omp` specifies how many omp threads to use for the parts of the code that can be executed in parallel. Since the classical schemes can usually be run efficiently with a single thread, parallelization is available only for the quantum schemes. Default `omp=1`

  * `qstls-fix` specifies the name of the binary file used to load the fixed component of the auxiliary density response function in the solution of the qSTLS scheme. If a file name is not specified, the fixed component of the auxiliary density response is computed from scratch and a significant increase in the computational cost for the solution of the qSTLS scheme can be expected. It is useful to note that state points with the same degeneracy parameter share the same fixed component of the auxiliary density response. Default `qslts-fix=NO_FILE` (no file is specified)

  * `qstls-guess` specifies the name of the binary file used as initial guess for the solution of the quantum schemes. If a file name is not specified, the iterative scheme is started by setting the auxiliary density response to zero (which is equivalent to set the dynamic local field correction to zero). Default `qslts-guess=NO_FILE` (no file is specified) 

  * `qstls-iet-fix` specifies the name of the folder containing the binary files used to load the fixed component of the auxiliary density response function in the solution of the qSTLS-IET scheme.  The files are assumed to have the name psi_fixed_theta\*_xx\**_.bin, where \* is the value of the quantum degeneracy parameter and \** is the value of the wave-vector. One file for each wave-vector is expected. If a folder name is not specified, the fixed component of the auxiliary density response is computed from scratch and a significant increase in the computational cost for the solution of the qSTLS-IET scheme can be expected. It is useful to note that state points with the same degeneracy parameter share the same fixed component of the auxiliary density response. Default `qslts-iet-fix=NO_FILE` (no folder is specified)

  * `--rs`  specifies the  quantum coupling parameter. Default `--rs=1.0`.
  
  * `--stls-guess` specifies the name of the binary file used as initial guess for the solution of the classical schemes. If a file name is not specified, the iterative scheme is started by setting the static local field correction to zero. Default `slts-guess=NO_FILE` (no file is specified) 
  
  * `--theory` specifies which scheme should be employed to compute the static structure factor. Default `--theory=STLS`. Accepted options are: 
     * `STLS` (classical STLS scheme), 
     * `STLS-HNC` (classical STLS-HNC scheme) 
     * `STLS-IET-IOI` (classical STLS-IET scheme with the [IOI bridge function](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.1051))
     * `STLS-IET-LCT` (classical STLS-IET scheme with the [LCT bridge function](https://arxiv.org/abs/2108.09574))
     * `QSTLS` (quantum STLS scheme). 
     * `QSTLS-IET-HNC` (quantum STLS-HNC scheme). 
     * `QSTLS-IET-IOI` (quantum STLS-IET scheme with the [IOI bridge function](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.1051))
     * `QSTLS-IET-LCT` (quantum STLS-IET scheme with the [LCT bridge function](https://arxiv.org/abs/2108.09574))
     
  * `--Theta` specifies the  quantum degeneracy parameter. Default `--Theta=1.0`
  
  * `--xmax` specifies the cutoff for wave-vector grid. Default `--xcut=20`
 
  ## Output 
  
  Regardless of the theory that is being solved, STLS produces the following output:
  
  * One text file with the static structure factor (ssf_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One text file with the static local field correction (slfc_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One text file with the static density response (sdr_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One text file with the ideal density response (idr_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One text file with the static structure factor within the Hartree-Fock approximation (ssfHF_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One text file with the interaction energy (uint_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One binary file which can be used as an initial guess for subsequent calculations via the options `--stls-guess` and `--qstls-guess` (restart_rs\*\_theta\*\*\_\*\*\*_.bin)

 For the quantum schemes, additional files are produced: 
 
  * One text file with the auxiliary density response (adr_rs\*\_theta\*\*\_\*\*\*_.dat)
  * One binary file which can be used to avoid recomputing the fixed component of the auxiliary density response in the qSTLS scheme via the option `--qstls-fix`  (fixed_rs\*\_theta\*\*\_\*\*\*.bin) **(only for the qSLTLS scheme)** 
  * One binary file per wave-vector which can be used to avoid recomputing the fixed component of the auxiliary density response in the qSTLS-IET scheme via the option `--qstls-iet-fix`  (psi_fixed_theta\*\*\_xx\*\*\*\*.bin, for N wave-vectors N files are written) **(only for the qSTLS-IET scheme)**
  
In the above \* is corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter, \*\*\* to the dielectric scheme that was solved and \*\*\*\* is the wave-vector value.

## Guess and restart

All the schemes solved by STLS are solved iteratively. In order to improve convergence, it is often beneficial to provide a suitable initial guess. This can be done via the options `stls-guess` and `qslts-guess`, that can be used to provide binary guess files which contain the initial guess (information on how to use these options is given in the section **Running**). The guess files for the classical schemes (STLS, STLS-HNC and STLS-IET) contain information regarding the static structure factor and the static local field correction, while the guess files for the quantum schemes (qSTLS, qSTLS-HNC and qSTLS-IET) contain information regarding the static structure factor and the auxiliary density response. Hence, guess files for classical and quantum scheme are not interchangeable. Obviously, the guess files can also be used as a restart point for all those calculations that were stopped before it was possible to achieve the desired level of accuracy (this can usually happen if the limit on the maximum number of iterations is hit before the error falls below the desider threshold).
