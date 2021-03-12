function out = PS_example23_alg
% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Diekmann, Vermiglio, Numerical bifurcation analysis of renewal
% equations via pseudospectral approximation, available at 
% https://arxiv.org/abs/2012.05364
  
%% PS_example23_alg.m
% MatCont system definition file of the PseudoSpectral Discretization of
% the nonlinear RE taken from Breda et al, EJQTDE 2016,
% x(t) = beta*int_a^amax x(t-a)*exp(x(t-a))da
% with the method of inversion of algebraic condition [EJQTDE 2016]
% using Chebyshev zeros (plus 0)
% The code uses the function polint.m, available from the Differentiation
% Matrix Suite by Weideman, Reddy, 2000

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= @userf; %user function to select specific parameter values
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,loggamma,abar,tau,M) 

    %% discretization of the unitary interval [-1,0]
    % construction of nodes and differentiation matrix
    p = pi*(2*(0:M-1)'+1)/(2*M);
    x=[1;sin(pi/2-p)]; % nodes with addition of 1 % either cos(p) or sin (pi/2-p) 
    X=repmat(x,1,M+1);
    dX=X-X';
    c=[2^(M-1)/M*prod(dX(1,2:end)); ((-1).^(0:M-1)').*dX(2:end,1)./sin(p)];
    D=(c*(1./c'))./(dX+(eye(M+1)));
    D=D-diag(sum(D')); % differentiation matrix

    % scaling
    Nodes = 0.5*tau*(x-1);
    DD = 2/tau*D;

    dMDM_RE = DD(2:end,:);

    %% SYSTEM DEFINITION *** specific to the equation ***

    % Options
    opt_fsolve = optimoptions('fsolve','Display','off','TolFun',1e-10); 
    %opt_ode = odeset('RelTol',1e-6,'AbsTol',1e-6); 

    % For quadrature formulas (efficiency)
    [QuadWeights,QuadNodes]=cheb_quad(50,-tau,-abar);

    f = @(x) x.*exp(-x);
    FM = @(x) 0.5*exp(loggamma)*QuadWeights*f(polint(Nodes,[x;state],QuadNodes));

    %% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

    if (abs(state(1)-FM(state(1)))==Inf || isnan(state(1)-FM(state(1))))
        KM = 0;
    else
        KM = fsolve(@(x) x-FM(x),state(1),opt_fsolve);
    end

    dydt= dMDM_RE*[KM;state];

end
 
 
% --------------------------------------------------------------------------
function Weq=init(M,xeq,yeq)
% INPUT: M is the discretization parameter
%        xeq,yeq are column vectors with the equilibrium states of the RE and
%        DDE respectively
% OUTPUT Weq is the initial vector for init_EP_EP
    
    Weq=[kron(ones(M,1),xeq); kron(ones(M+1,1),yeq)];
    
end

function out=userf(time,state,loggamma,abar,aux,tau,M) 
% Userfunction to select specific value of parameter
    
    out = (loggamma-1).*(loggamma-2).*(loggamma-3).*(loggamma-4);
    
end

% ------------
%%% AUXILIARY FUNCTIONS
function [w,x]=cheb_quad(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% see Trefethen 2000

    p=pi*(0:N)'/N;
    x=((a-b)*cos(p)+b+a)/2;

    % Quadrature weights
    w=zeros(1,N+1);
    ii=2:N;
    v=ones(N-1,1);
    if mod(N,2)==0
        w(1)=1/(N^2-1);
        w(N+1)=w(1);
        for k=1:N/2-1
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
        v=v-cos(N*p(ii))/(N^2-1);
    else
        w(1)=1/N^2;
        w(N+1)=w(1);
        for k=1:(N-1)/2
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
    end
    w(ii)=2*v/N;
    w=w*abs(b-a)/2;

end


