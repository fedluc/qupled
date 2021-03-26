# STLS

## Introduction

stls computes the static structure factor and thermodynamic properties of quantum one component plasmas via theoretical approaches based on the dielectric formalism. The theoretical approaches which can be solved with stls include:

* The classical STLS approach as discussed by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278)
* The classical STLS-HNC approach as discussed by [Tanaka](https://aip.scitation.org/doi/full/10.1063/1.4969071)
* The quantum STLS appraoach as discussed by [Schweng](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037)

## Limitations

stls can only be employed to compute the static structure factor of finite-temperature systems. Calculations of ground state properties will be implemented in the future.

## Compiling

The code can be compiled with the [make file](src/Makefile) provided in the source directory. Please note that in order to correctly compile the program it is necessary 
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

stls accepts the following command line options (the same information can also be retrieved by running stls with the option `--help`) :

  * `-d` or `--dx` specifies the  resolution for wave-vector grid. Default `-d 0.01` or `--dx=0.01`
  * `-e` or `--errIter` specifies the minimum error for convergence in the iterative procedure used to compute the static structure factor.  Default `-e 1e-5` or `--errIter=1e-5`
  * `-g` or `--mg=-10,10` specifies the  initial guess for chemical potential
  * `-i` or `--iter=1000`  specifies the maximum number of iterations to employ in
  * `-l` or `--nl=128`   specifies the   number of Matsubara frequencies
  * `-m` or `--mix=0.1` specifies the             mixing parameter for iterative solution
  * `-p` or `--phi=NO_FILE` specifies the        load density response from PHI_FILE
  * `-r` or `--rs=1.0`    specifies the           quantum coupling parameter
  * `-s` or `--sol=STLS`   specifies the          theory to be solved
  * `-t` or `--Theta=1.0` specifies the           quantum degeneracy parameter
  * `-u` or `--uex=NO_FILE` specifies the         compute internal energy from data in SSF_FILE
  * `-x` or `--xcut=20.48` specifies the          cutoff for wave-vector grid
 

  
  
  ## Output 
  
  The output of the code consists of:
  
  * One text file with the static structure factor
  * One text file with the static local field correction
  * One binary file with the density response. Since the density
        response depends only on Theta, this file can be stored and
        provided in input for subsequent solutions of the STLS approach
        with the same Theta (see option `-p`). It should be noted
        that, if the option -p is used, the values of the quantum
        degeneracy parameter, of the grid resolution and of the grid
        cutoff specified in input will be overwritten by the values
        contained in the density response file provided in input
-->
