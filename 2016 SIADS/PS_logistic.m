function out = PS_logistic
% Matcont system definition file of the PseudoSpectral Discretization of
% the delayed logistic equation
% y'(t) = r*y(t)*[1-y(t-1)]

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= []; %@delay; % function computing the delay, given the state vector
out{11}= []; %@userf1;
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,r,aux,tau,M) 
% OUTPUT: evaluation of the rhs of the approximating ODE system
% INPUT: (notation as in Breda et al., SIADS 2016)
% state=(UM,yM,VM) vector of variables, dimension M*d1+(M+1)*d2
% M: discretization index
% aux: auxiliary parameter (e.g. for BP continuation)
% Notice that: parameters MUST be listed separately

%% PHASE 0:
% discretization of the unitary interval [-1,0]

[UnitQuadweights,UnitNodes,UnitDD,BaryWeights]=cheb(M,-1,0);

%% PHASE 1: SYSTEM DEFINITION *** to be completed by the user ***
% It must include definition of:
% d1, d2: dimension of the RE and DDE, respectively
% tau_max: maximal fixed delay for the discretization interval [-tau_max,0]
% FM, GM: rhs of the RE and DDE, respectively. N.B.: they must be defined
% as function of x=x(t) (state of the RE in theta_0=0)
% to interpolate the state [yM,VM] in a point theta \neq 0,-tau_max,
% use barycentric interpolation:
% interpoly(theta,tau_max*UnitNodes,[yM;VM],BaryWeights)

% Options
%opt_fsolve = optimoptions('fsolve','Display','off','TolFun',1e-6); 

% Parameters and functions
d1=0; % dimension of the RE
d2=1; % dimension of the DDE
tau_max=tau;

%UM=state(1:d1*M);
yM=state((d1*M+1):(d1*M+d2)); VM=state((d1*M+d2+1):end);
FM = @(x) [];
GM = @(x) r*yM.*(1-VM(end));

%% FINAL PHASE
% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

KM=[]; 
if d1>0
    KM = fsolve(@(x) x-FM(x),UM(1:d1),opt_fsolve);
end
%dMDM_RE=kron(UnitDD(2:end,:),eye(d1));
dMDM_DDE=kron(UnitDD(2:end,:),eye(d2));

dydt= [ %(1/tau_max*dMDM_RE)*[KM;UM];
        GM(KM);(1/tau_max*dMDM_DDE)*[yM;VM]];
  
end
 
 
% --------------------------------------------------------------------------
function state_eq=init(M,xeq,yeq)
% INPUT: M is the discretization parameter
%        xeq,yeq are column vectors with the equilibrium states of the RE and
%        DDE respectively
% OUTPUT state_eq is the initial vector for init_EP_EP
    
    state_eq=[kron(ones(M,1),xeq); kron(ones(M+1,1),yeq)];
    
end

% ------------
%%% AUXILIARY FUNCTIONS

function [w,x,D,q]=cheb(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=b, x_N=a),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
% see Trefethen

if N==0
    x=1;
    D=0;
    return
end
p=pi*(0:N)'/N;
x=((b-a)*cos(p)+b+a)/2;

c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=x(:,ones(1,N+1)); %X=repmat(x,1,N+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D,2)); %D=D-diag(sum(D'));

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

% Barycentric weights
q=1./prod(dX'+eye(N+1)); %q=1./prod(dX'+eye(N+1)); % row vector of the barycentric weights
end