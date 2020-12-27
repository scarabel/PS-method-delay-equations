READ ME

This repository contains some codes that are useful to reproduce the results in the papers

[SIADS2016] Breda D., Diekmann O., Gyllenberg M., Scarabel F., Vermiglio R. (2016). Pseudospectral discretization of nonlinear delay equations: new prospects for numerical bifurcation analysis, SIAM Journal on applied dynamical systems, 15(1), 1–23.

Each example consists of two files:
PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
MC_example: script for the Matcont continuation of the system defined in "PS_example".

The codes are tested on Matlab 2014a and Matcont version cl_matcont5p4, 
and on Matlab 2017a and Matcont version matcont6p6.

## Instruction for the PS_example file
"PS_example" is the system definition file required by Matcont. 
It contains the definition of the right-hand side of the PseudoSpectral discretization of the system under study.
To perform the Matcont continuation, the system definition file ``PS_example'' must be copied into the subfolder "systems" of the Matcont folder.


   function dydt = fun_eval(time,state,p1,...,aux,tau,M) 
   
   OUTPUT: evaluation of the rhs of the approximating ODE system
   INPUT: (notation as in Breda et al., SIADS 2016)
   state=(UM,yM,VM) vector of variables, dimension M*d1+(M+1)*d2
   M: discretization index
   aux: auxiliary parameter (e.g. for BP continuation)
   p1,...: model parameters (specific of the model)
   Notice that: model parameters MUST be listed separately
   
   The user should modify only the section
   PHASE 1: SYSTEM DEFINITION *** to be completed by the user ***
   specifying the desired options and the parameters and functions of the model
   
   PHASE 1 must include the definition of:
   d1, d2: dimension of the RE and DDE, respectively
   tau_max: maximal fixed delay for the discretization interval [-tau_max,0]
   FM, GM: rhs of the RE and DDE, respectively. N.B.: they must be defined as function of x=x(t) (state of the RE in theta_0=0)
   Note that, to interpolate a vector W (dimension M+1) at a point theta \neq 0,-tau_max, we use barycentric interpolation:
       interpoly(theta,tau_max*UnitNodes,W,BaryWeights)
   
   Options:
   ``opt_fsolve'' is used to invert the algebraic equation coming from RE
   ``opt_ode'' possible options for solving ode systems %opt_ode = odeset('RelTol',1e-6,'AbsTol',1e-6); 

The user may also modify the section FINAL PHASE in order to treat RE or the DDE only, for speed up reason.
However, this is not necessary.


## Instructions for the MC_example file
"MC_example" is the script containing the instruction for the Matcont continuation.
Continuation and parameter options may be varied in the script.
