function out = PS_system
% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O. (2020)
% Pseudospectral approximation of Hopf bifurcation for delay differential
% equations, SIAM Journal on Applied Dynamical Systems.
%
%% PS_blowflies.m
% MatCont system definition file of the PseudoSpectral Discretization of
% the neural system
% w'(t) = 1-0.5*k*w(t)*w(t-1)*q(t-1)
% q'(t) = w(t)-c
% Using Chebyshev zeros (plus 0 and -1)
% The code uses the code poldif.m from the Differentiation Matrix Suite
% (Weideman, Reddy, 2000)

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= @userf; %@delay; % function computing the delay, given the state vector
out{11}= []; %@userf1;
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,k,c,M) 

%% Discretization of the unitary interval [-1,0]

tau = 1;
% define Chebyshev zeros
angles = pi*(2*[1:M-1]'-1)/(2*(M-1));
Nodes = [0;0.5*cos(angles)-0.5;-tau];
% compute differentiation matrix using Weideman and Reddy
DD = poldif(Nodes,1);
D = kron(DD(2:end,:),eye(2));

%% SYSTEM DEFINITION *** to be completed by the user ***

% Parameters and functions

W = state(1:2:end); % discrete state vector for variable w(t)
Q = state(2:2:end); % discrete state vector for variable q(t)

dydt= [ 
    1 - 0.5*k*W(1)*W(end)*Q(end);
    W(1) - c;
    D*state];

end
 
 
% --------------------------------------------------------------------------
function Weq=init(M,yeq)
% INPUT: M is the discretization parameter
%        yeq is the column vector of the equilibrium state
% OUTPUT Weq is the initial vector for init_EP_EP
    
    Weq=kron(ones(M+1,1),yeq);
    
end

function y = userf(time,state,k,c,M)
    y=(k-1.5);
end
