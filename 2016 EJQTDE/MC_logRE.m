% MC_logRE
% Matcont continuation of the 
% nonlinear RE with quadratic nonlinearity defined in PS_logRE

clear;
clearvars -global cds

% Discretization parameters
M=10;

% Initial parameter values
d1=1; d2=0; % dimension of RE and DDE, respectively
tau_max=3;
abar=1;
gamma=1.2;
epsi=1;
aux=1;
par=[gamma,abar,epsi,aux,tau_max,M]';

% Approximated equilibrium corresponding to par
xeq=1-1/gamma;
yeq=[];

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
%ap3=3;
TOL=1e-6;

%% Continuation process

MM=d1*M+d2*(M+1); % dimension of the approximating ODE system
handles=feval(@PS_logRE);
opt=contset;
global cds;

%% Equilibrium continuation from initial point [xeq;yeq]

display('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
%opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',50);

Weq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_logRE,Weq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< 5)
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
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_logRE,H,parH,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< 4)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
figure(1)
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
xlabel('parameter','interpreter','latex');
%ylabel('max/min','interpreter','latex')

%% Detection of singular points
% PD, Period Doubling (on a branch of periodic solutions)
for ii=1:size(slc)
    if strcmp(slc(ii).label,'PD ')==1
        PD_index=slc(ii).index;
        sPD=slc(ii);
        break;
    end
end
par(ap1)=xlc(end,PD_index);
PD=xlc(1:MM,PD_index);

xlcPD=xlc; vlcPD=vlc; slcPD=slc; hlcPD=hlc; flcPD=flc;
parPD=par;

%% Limit cycle continuation from PD
% xlcPD = output of the continuation where PD was detected
% parPD = parameter vector at PD
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from PD1');

% set options
opt=contset(opt,'MaxStepSize',0.001);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_PD_LC(@PS_MackeyGlass,xlcPD,sPD,ntst,ncol,1e-6);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<4) && xlc(end,end)< 4)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
% figure
upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
xlabel('parameter','interpreter','latex');
ylabel('max/min','interpreter','latex')
% or equivalently my own function cpl_per:
% cpl_per(xlc,slc,flc,ntst,ncol,MM,M+1,size(xlc,2));

% Plot profile of orbit corresponding to index j
% j='end';
% figure
% mesh= xlc(end-1,j)*flc(1:ntst+1,j);
% profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),j);
% plot(mesh,profile);

%% Detection of singular points
% PD, Period Doubling (on a branch of periodic solutions)
for ii=1:size(slc)
    if strcmp(slc(ii).label,'PD ')==1
        PD_index=slc(ii).index;
        sPD=slc(ii);
        break;
    end
end
par(ap1)=xlc(end,PD_index);
PD=xlc(1:MM,PD_index);

xlcPD=xlc; vlcPD=vlc; slcPD=slc; hlcPD=hlc; flcPD=flc;
parPD=par;

%% Limit cycle continuation from PD
% xlcPD = output of the continuation where PD was detected
% parPD = parameter vector at PD
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from PD2');

% set options
opt=contset(opt,'Adapt',0);
TOL=1e-3;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_PD_LC(@PS_MackeyGlass,xlcPD,sPD,ntst,ncol,1e-6);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< 4)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
% figure
upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
xlabel('parameter','interpreter','latex');
ylabel('max/min','interpreter','latex')
% or equivalently my own function cpl_per:
% cpl_per(xlc,slc,flc,ntst,ncol,MM,M+1,size(xlc,2));

% Plot profile of orbit corresponding to index j
% j='end';
% figure
% mesh= xlc(end-1,j)*flc(1:ntst+1,j);
% profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),j);
% plot(mesh,profile);

%% Detection of singular points
% PD, Period Doubling (on a branch of periodic solutions)
for ii=1:size(slc)
    if strcmp(slc(ii).label,'PD ')==1
        PD_index=slc(ii).index;
        sPD=slc(ii);
        break;
    end
end
par(ap1)=xlc(end,PD_index);
PD=xlc(1:MM,PD_index);

xlcPD=xlc; vlcPD=vlc; slcPD=slc; hlcPD=hlc; flcPD=flc;
parPD=par;
