# 2021 SIADS Hopf

This folder contains some codes that are useful to reproduce the results in the paper

De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O.
Pseudospectral approximation of Hopf bifurcation for delay differential equations, 
SIAM Journal on Applied Dynamical Systems. Preprint available at https://arxiv.org/abs/2006.13810 

The codes are used to perform the numerical bifurcation of the blowflies equation and the system for neural network.
The numerical bifurcation is performed using the following Matlab packages:
- dde-biftool, availablet at: https://sourceforge.net/projects/ddebiftool/
- MatCont, available at: https://sourceforge.net/projects/matcont/

The differentiation matrices are computed using the function poldif.m contained in the Differentiation Matrix Suite
from Weideman, Reddy 2000: http://appliedmaths.sun.ac.za/~weideman/research/differ.html

For the Matcont continuation, each example consists of two files:
1) PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
To perform the Matcont continuation, the system definition file "PS_example" must be copied into the subfolder "systems" of the Matcont folder.

2) MC_example: script for the Matcont continuation of the system defined in "PS_example".

The codes are tested on MATLAB 2019b and MatCont version MatCont7p1.

