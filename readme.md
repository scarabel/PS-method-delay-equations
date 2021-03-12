# PS-method-delay-equations
 Pseudospectral discretization of delay equations

This repository contains the MATLAB codes for the pseudospectral discretization of delay equations (delay differential and renewal equations), and the numerical bifurcation analysis using the package MatCont for MATLAB.

The main references are
- for delay equations (DDE and renewal equations):
[SIADS2016] Breda D, Diekmann O, Gyllenberg M, Scarabel F, Vermiglio R (2016). Pseudospectral discretization of nonlinear delay equations: new prospects for numerical bifurcation analysis, SIAM Journal on applied dynamical systems, 15(1), 1–23. https://doi.org/10.1137/15M1040931 
- for renewal equations (original method with inversion of the algebraic equation)
[EJQTDE2016] Breda D, Diekmann O, Liessi D, Scarabel F (2016). Numerical bifurcation analysis of a class of nonlinear renewal equations, Electronic Journal of Qualitative Theory of Differential Equations, 65, 1–24. https://doi.org/10.14232/ejqtde.2016.1.65 

Each example consists of two files:
PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the MatCont continuation.
MC_example: script for the MatCont continuation of the system defined in "PS_example".

The file Matcont_routines provides a collection of "building blocks" to help the user in coding their own MatCont continuation.
This does not cover all the possibilities of MatCont, but can be used to extract and adapt from time to time only the sections of the code suitable for performing 
the continuation of the specific system of interest.

The codes are released under the MIT license (see file LICENSE.txt for details).
