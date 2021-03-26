# STLS

## Introduction

stls computes the static structure factor and thermodynamic properties of quantum one component plasmas via theoretical approaches based on the dielectric formalism. The theoretical approaches which can be solved with stls include:

* The classical STLS approach as discussed by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278)
* The classical STLS-HNC approach as discussed by [Tanaka](https://aip.scitation.org/doi/full/10.1063/1.4969071)
* The quantum STLS appraoach as discussed by [Schweng](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037)

## Limitations

stls can only be employed to compute the static structure factor of finite-temperature systems. Calculations of ground state properties will be implemented in the future.

## Compiling

The code can be compiled with gcc and with the [make file](src/Makefile) provided in the source directory. Please note that in order to correctly compile the program it is necessary 
that the [GNU scientific library](https://www.gnu.org/software/gsl/) is installed and that the path to the library is included in the make file.

## Running 

Given a state point defined via the quantum degeneracy parameter (Theta) and via the quantum coupling parameter (r<sub>s</sub>), stls computes the static structure factor with the following procedure:

* The chemical potential is determined from the normalization condition for the Fermi-Dirac distribution function which is solved with a bisection method
* The normalized Lindhard density response is computed for various values of the wave-vector and of the Matsubara frequency.
* For the STLS and STLS-HNC approaches, the static structure factor together and the static local field correction are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||G<sub>i</sub>(x) - A<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here G(x) is the static local field correction and epsilon is a tolerance specified in input.
* For the quantum STLS, the static structure factor and the dynamic field correctionis are computed via an iterative solution which employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
||S<sub>i</sub>(x) - S<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here S(x) is the static structure factor and epsilon is a tolerance specified in input
* The internal energy is calculated by integrating the static structure factor

All the integrals which appear in the calculation of the static structure factor and of the internal energy are computed with a mid-point Riemann sum.

The following command line options can be employed to control the calculations performed by stls (the same information can also be retrieved by running stls with the option `--help`) :

  * `-d` or `--dx` specifies the  resolution for wave-vector grid. Default `-d 0.01` or `--dx=0.01`
  * `-e` or `--errIter` specifies the minimum error for convergence in the iterative procedure used to compute the static structure factor.  Default `-e 1e-5` or `--errIter=1e-5`
  * `-g` or `--mg` specifies the initial guess for the bisection procedure used to compute the chemical potential.  Default `-g -10,10` or `--mg=-10,10`
  * `-i` or `--iter`  specifies the maximum number of iterations to employ in the iterative procedure used to compute the static structure factor.  Default `-i 1000` or `--iter=1000`
  * `-l` or `--nl=128`  specifies the number of Matsubara frequencies.  Default `-l 128` or `--nl=128`.
  * `-m` or `--mix=0.1` specifies the mixing parameter for iterative solution.  Default `-m 0.1` or `--mix=0.1`.
  * `-p` or `--phi=FILE` specifies to load the normalized Lindhard density response from a file instead of computing it. The file with the density response must be a binary file produced by a previous run of stls (see the output section). Default `-p NO_FILE` or `--phi=NO_FILE` (no file is loaded and the normalized Lindhard density is computed from scratch).
  * `-r` or `--rs=1.0`    specifies the  quantum coupling parameter. Default `-r 1.0` or `--rs=1.0`.
  * `-s` or `--sol=STLS`   specifies which approach should be employed to compute the static structure factor. Accepted options are `STLS` (classical STLS approach), `STLS-HNC` (classical STLS-HNC approach) and `QSTLS` (quantum STLS approach). Default `-s STLS` or `--sol=STLS`.
  * `-t` or `--Theta=1.0` specifies the  quantum degeneracy parameter. Default `-t 1.0` or `--Theta=1.0`.
  * `-x` or `--xcut=20.48` specifies the cutoff for wave-vector grid. Default `-x 20.48` or `--xcut=20.48`.
 
  ## Output 
  
  stls produces the following output:
  
  * One text file with the static structure factor
  * One text file with the static local field correction 
  * One text file with the dynamic local field correction (only for the quantum STLS approach)
  * One binary file with the density response. Since the density response depends only on Theta, this file can be stored an provided in input for subsequent solutions of the STLS approach with the same Theta (see option `-p`). It should be noted that, if the option -p is used, the values of the quantum degeneracy parameter, of the grid resolution and of the grid cutoff specified in input will be overwritten by the values contained in the density response file provided in input

## GPU 

An implementation of the quantum STLS approach which can be run on graphics cards (GPU) is available in the folder [MOD_GPU](src/MOD_GPU). For the default parameters employed in stls, the GPU implementation is approximately two order of magnitude faster than the corresponding CPU version. However, all the calculations in the GPU implementation are performed in single precision. The GPU implementation can be compiled with a [dedicated make file](src/Makefile) and requires the nvcc compiler.
