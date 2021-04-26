# STLS

## Introduction

stls computes the static structure factor and thermodynamic properties of quantum one component plasmas via theoretical approaches based on the dielectric formalism. The theoretical approaches which can be solved with stls include:

* The classical STLS approach as discussed by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278)
* The classical STLS-HNC approach as discussed by [Tanaka](https://aip.scitation.org/doi/full/10.1063/1.4969071)
* The classical STLS-IET approach which is obtained by complementing the classical STLS-HNC approach with the OCP bridge function parameterization of [Ichimaru and collaborators](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.46.1051)

## Limitations

stls can only be employed to compute the static structure factor of finite-temperature systems. Calculations of ground state properties will be implemented in the future.

## Compiling

The code can be compiled with gcc and with the [make file](Makefile) provided in the source directory. Please note that in order to correctly compile the program it is necessary that the following libraries are installed

* [GNU scientific library](https://www.gnu.org/software/gsl/). This library must be explicitly included by calling `make GSL="PATH"`, where `PATH` is the path to the folder containing the header files of the library. Alternatively, it is possible to change the default value of the GSL variable in the Makefile and then to simply compile via `make`.
* [OpenMP library](https://en.wikipedia.org/wiki/OpenMP). In most cases it is not necessary to explicitly include the path to this library. However, if this should be necessary, it should be possible to do so by modifying the `INCLUDE` variables which appear in the Makefiles in the [cquad](cquad) and [riemann](riemann) folders.

## Running 

Given a state point defined via the quantum degeneracy parameter (Theta) and via the quantum coupling parameter (r<sub>s</sub>), stls computes the static structure factor with the following procedure:

* The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method
* The normalized Lindhard density response is computed for various values of the wave-vector and of the Matsubara frequency.
* The static structure factor and the static local field correction are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||G<sub>i</sub>(x) - G<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here G(x) is the static local field correction and epsilon is a tolerance specified in input.
* The internal energy is calculated by integrating the static structure factor

All the integrals which appear in the calculation of the static structure factor and of the internal energy are computed with a mid-point Riemann sum.

The following command line options can be employed to control the calculations performed by stls (the same information can also be retrieved by running stls with the option `--help`) :

  * `-d` or `--dx` specifies the  resolution for wave-vector grid. Default `-d 0.01` or `--dx=0.01`
  * `-e` or `--errIter` specifies the minimum error for convergence in the iterative procedure used to compute the static structure factor.  Default `-e 1e-5` or `--errIter=1e-5`
  * `-f` or `--sg` specifies to load the initial guess from a file. The file with the initial guess must be a binary file produced by a previous run of stls (see the output section). Default `-f NO_FILE` or `--sg=NO_FILE` (no file is loaded and the initial guess is assumed to be given by the random phase approximation solution, i.e. zero static local field correction).
  * `-g` or `--mg` specifies the initial guess for the bisection procedure used to compute the chemical potential.  Default `-g -10,10` or `--mg=-10,10`
  * `-i` or `--iter`  specifies the maximum number of iterations to employ in the iterative procedure used to compute the static structure factor.  Default `-i 1000` or `--iter=1000`
  * `-l` or `--nl=128`  specifies the number of Matsubara frequencies.  Default `-l 128` or `--nl=128`.
  * `-m` or `--mix=0.1` specifies the mixing parameter for iterative solution.  Default `-m 0.1` or `--mix=0.1`.
  * `-o` or `--omp=1` specifies how many omp threads to use for the parts of the code that can be executed in parallel. For single core calculations, 16 threads is usually a reasonable choice. Default `-o 1` or `--omp=1`.
  * `-r` or `--rs=1.0`    specifies the  quantum coupling parameter. Default `-r 1.0` or `--rs=1.0`.
  * `-s` or `--sol=STLS`   specifies which approach should be employed to compute the static structure factor. Accepted options are `STLS` (classical STLS approach), `STLS-HNC` (classical STLS-HNC approach) and `STLS-IET` (classical STLS-IET approach). Default `-s STLS` or `--sol=STLS`.
  * `-t` or `--Theta=1.0` specifies the  quantum degeneracy parameter. Default `-t 1.0` or `--Theta=1.0`.
  * `-x` or `--xcut=20.48` specifies the cutoff for wave-vector grid. Default `-x 20.48` or `--xcut=20.48`.
 
  ## Output 
  
  stls produces the following output:
  
  * One text file with the static structure factor (ssf_rs\*_theta\*\*\_\*\*\*.dat)
  * One text file with the static local field correction (slfc_rs\*_theta\*\*\_\*\*\*.dat)
  * One text file with the static density response (sdr_rs\*_theta\*\*\_\*\*\*.dat)
  * One text file with the normalized ideal Lindhard density response (idr_rs\*_theta\*\*\_\*\*\*.dat)
  * One text file with the static structure factor within the Hartree-Fock approximation (ssfHF_rs\*_theta\*\*\_\*\*\*.dat)
  * One binary file which can be used as an initial guess for subsequent calculations via the option `-f` (restart_rs\*_theta\*\*\_\*\*\*.bin)

In the above \* is corresponds to the value of the quantum coupling parameter, \*\* to the value of the quantum degeneracy parameter and \*\*\* to the dielectric approach that was solved. 

## Quantum STLS 

The folder [QSTLS](QSTLS) containst an implementation of the quantum STLS appraoach as discussed by [Schweng](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037). Due to recent updates in the stls code, the QSTLS implementation is temporarily not compatible with the other parts of the stls code and it is therefore confined in a dedicated folder. The QSTLS will be incorporated in the stls code in the future (updated: April 21, 2021)

An implementation of the quantum STLS approach which can be run on graphics cards (GPU) is available in the folder [MOD_GPU](MOD_GPU). For the default parameters employed in stls, the GPU implementation is approximately two order of magnitude faster than the corresponding CPU version. However, all the calculations in the GPU implementation are performed in single precision. The GPU implementation can be compiled with a [dedicated make file](MOD_GPU/Makefile) and requires the nvcc compiler.
