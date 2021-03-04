# STLS

stls solves the classical STLS approach as defined in by [Tanaka and Ichimaru](https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278). The state 
  point of interest is defined via the quantum degeneracy parameter (Theta)
  and via the quantum coupling parameters (r<sub>s</sub>). The equation for the 
  chemical potential is solved via bisection method for which two 
  initial guesses must be provided in input via the option `-g`.
  The STLS approach is solved iteratively on a wave-vector grid 
  extending from 0 to agiven cutoff. The grid resolution and cutoff
  are specified in input together with number of Matsubara frequencies
   necessary for the calculation of the static structure factor. The 
  iterative solution employs [mixing](https://aip.scitation.org/doi/abs/10.1063/1.1682399]) and is assumed converged if the condition 
  ||G<sub>i</sub>(x) - G<sub>i-1</sub>(x)|| < epsilon is satisfied between two successive iterations. Here G(x) is the 
  static local field correction, epsilon is a tolerance specified in
  input. The output of the code consists of:
  
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
