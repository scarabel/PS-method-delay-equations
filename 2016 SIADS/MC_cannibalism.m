% MC_cannibalism
% Matcont continuation of the 
% cannibalism model RE defined in PS_cannibalism

clear;
clearvars -global cds

% Discretization parameters
M=20;

% Initial parameter values
d1=1; d2=0; % dimension of RE and DDE, respectively
tau_max=4;
abar=3;
beta=0;
aux=1;
par=[beta,abar,aux,tau_max,M]';

% Approximated equilibrium corresponding to par
xeq=0;
yeq=[];

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
%ap3=3;
TOL=1e-10;
TestTOL=1e-6;

%% Continuation process

MM=d1*M+d2*(M+1); % dimension of the approximating ODE system
handles=feval(@PS_cannibalism);
opt=contset;
global cds;

%% Equilibrium continuation from initial point [xeq;yeq]

display('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TestTOL);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',50);

state_eq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_cannibalism,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< 10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
hold on;

%% Detection of singular points
% xe,ve,se,he,fe = output of previous continuation
% par = current parameter vector
% ap1 = index of continuation parameter in vector par

% BP, branching point
for ii=1:length(se)
    if strcmp(se(ii).label,'BP')==1
        BP_index=se(ii).index;
        sBP=se(ii);
        break;
    end
end
par(ap1)=xe(end,BP_index);
BP=xe(1:MM,BP_index);

xeBP=xe; veBP=ve; seBP=se; heBP=he; feBP=fe;
parBP=par;

%% Equilibrium continuation from BP
% BP = vector of variables at BP
% parBP = parameter vector at BP
% sBP = information about BP from previous continuation
display('Starting equilibrium continuation from BP');

% set options
%opt=contset(opt,'MaxNumPoints',50);

[x0,v0]=init_BP_EP(@PS_cannibalism,BP,parBP,sBP,0.001);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
end

% Plot
figure(1)
cpl(xe,ve,se,[MM+1 1]);

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
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',0);
opt=contset(opt,'Multipliers',1);

ntst=40; % number of intervals
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_cannibalism,H,parH,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< 8)
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
j=257;
figure
mesh= xlc(end-1,j)*flc(1:ntst+1,j);
profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),j);
plot(mesh,profile,'r');