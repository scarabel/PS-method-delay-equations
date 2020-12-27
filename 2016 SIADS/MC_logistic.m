% MC_logistic
% Matcont continuation of the 
% delayed logistic equation defined in PS_logistic

clear;
clearvars -global

% Discretization parameters
M=10;

% Initial parameter values
d1=0; d2=1; % dimension of RE and DDE, respectively
tau_max=1;
r=0.1;
aux=1; % auxiliary parameter for codim-2 bifurcations
par=[r,aux,tau_max,M]';

% Approximated equilibrium corresponding to par
xeq=[];
yeq=1;

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
%ap3=3;
TOL=1e-10;

%% Continuation process

MM=d1*M+d2*(M+1); % dimension of the approximating ODE system
handles=feval(@PS_logistic);
opt=contset;
global cds;

%% Equilibrium continuation from initial point [xeq;yeq]

display('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);

state_eq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_logistic,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< 2)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
hold on;

%% Detection of singular points
% H, Hopf point

for ii=1:size(se)
    if strcmp(se(ii).label,'H ')==1
        H_index=se(ii).index;
        break;
    end
end
par(ap1)=xe(end,H_index);
H=xe(1:MM,H_index);

xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
parH=par;

%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

% set options
opt=contset(opt,'Singularities',0);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_logistic,H,parH,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< 2)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
figure(1)
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
xlabel('parameter','interpreter','latex');
%ylabel('max/min','interpreter','latex')

% Plot profile of orbit corresponding to index j
j=1;
figure
mesh= xlc(end-1,j)*flc(1:ntst+1,j);
profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),j);
plot(mesh,profile);