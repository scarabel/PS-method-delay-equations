# 2021 Renewal submitted

This folder contains some codes that are useful to reproduce the results in the paper
Scarabel, Diekmann, Vermiglio, Numerical bifurcation analysis of renewal equations via pseudospectral approximation, available at: https://arxiv.org/abs/2012.05364

The codes are used to perform the numerical bifurcation analysis of the examples in the paper.
The numerical bifurcation is performed using the Matlab package MatCont, available at: https://sourceforge.net/projects/matcont/

The interpoaltion is performed using polint.m contained in the Differentiation Matrix Suite
from Weideman, Reddy 2000: http://appliedmaths.sun.ac.za/~weideman/research/differ.html

For the Matcont continuation, each example consists of two files:
1) PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
2) MC_example: script for the MatCont continuation of the system defined in "PS_example".

The codes are tested on MATLAB 2020b and MatCont version MatCont7p1.
