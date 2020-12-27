function out = PS_blowflies
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
% Nicholson's blowflies equation
% y'(t) = -mu*y(t)+beta*y(t-1)*exp(-y(t-1))
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
function dydt = fun_eval(time,state,log_beta,mu,M) 

%% Discretization of the unitary interval [-1,0]

tau = 1;
% define Chebyshev zeros
angles = pi*(2*[1:M-1]'-1)/(2*(M-1));
Nodes = [0;0.5*cos(angles)-0.5;-tau];
% compute differentiation matrix using Weideman and Reddy, Differentiation
% Matrix Suite (2000)
DD = poldif(Nodes,1);
D = DD(2:end,:);

%% SYSTEM DEFINITION *** to be completed by the user ***

% Parameters and functions

dydt = [ -mu*state(1) + exp(log_beta).*state(end).*exp(-state(end));
        D*state];

end
 
 
% --------------------------------------------------------------------------
function Weq=init(M,yeq)
% INPUT: M is the discretization parameter
%        yeq is the scalar initial equilibrium guess 
% OUTPUT Weq is the initial vector for init_EP_EP
    
    Weq=yeq*ones(M+1,1);
    
end

function y = userf(time,state,log_beta,mu,M)
% to plot eigenvalues and multipliers
%    y=(log_beta-log(15)).*(log_beta-log(45)); % for equilibria (mu=3)
    y=log_beta-log(15*7); % for periodic solutions (mu=5)
end
