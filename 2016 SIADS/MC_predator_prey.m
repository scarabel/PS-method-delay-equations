% MC_logistic
% Matcont continuation of the 
% delayed logistic equation defined in PS_logistic

clear;
clearvars -global cds

% Discretization parameters
M=10;

% Initial parameter values
d1=0; d2=2; % dimension of RE and DDE, respectively
death=0.2;
juv_death=0.8;
conv=1;
pred=5;
aux=1; % auxiliary parameter for codim-2 bifurcations
tau=5;

par=[death,juv_death,conv,pred,aux,tau,M]';

% Approximated equilibrium corresponding to par
xeq=[];
yeq=[1;0];

% Continuation parameters
ap1=6; % index of the continuation parameter in the vector par
ap2=2;
ap3=5;
TOL=1e-6;

%% Continuation process

MM=d1*M+d2*(M+1); % dimension of the approximating ODE system
handles=feval(@PS_predator_prey);
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
opt=contset(opt,'Backward',1);

Weq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_predator_prey,Weq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)>0)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
axis([0 5 0 1.2])
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
opt=contset(opt,'Backward',1);

[x0,v0]=init_BP_EP(@PS_predator_prey,BP,parBP,sBP,0.001);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<5) && xe(end,end)>0)
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
opt=contset(opt,'Singularities',0);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_predator_prey,H,parH,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && min(xlc(end,:))>0.3)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
figure(1)
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
xlabel('parameter','interpreter','latex');
%ylabel('max/min','interpreter','latex')

%% BP continuation in two parameters
% BP = vector of variables at BP
% parBP = parameter vector at BP
% ap1,ap2,ap3 = index of continuation parameter in par (often ap3 is
% auxiliary parameter in degenerate biological models)
display('Starting BP continuation');

% set options
opt=contset(opt,'InitStepSize',0.001);
opt=contset(opt,'Backward',0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',0);

[x0,v0]=init_BP_BP(@PS_predator_prey,BP,parBP,[ap1 ap2 ap3],2);
[xbp,vbp,sbp,hbp,fbp]=cont(@branchpoint,x0,[],opt); xbp(MM+1,end)
while (xbp(MM+1,end)<5 && xbp(MM+1,end)>1)
    [xbp,vbp,sbp,hbp,fbp]=cont(xbp,vbp,sbp,hbp,fbp,cds); xbp(MM+1,end)
end

% Plot
figure(2); clf;
cpl(xbp,vbp,sbp,[MM+1 MM+2]);
hold on
% or equivalently:
% plot(xbp(MM+1,:),xbp(MM+2,:))

%% BP continuation in two parameters
% BP = vector of variables at BP
% parBP = parameter vector at BP
% ap1,ap2,ap3 = index of continuation parameter in par (often ap3 is
% auxiliary parameter in degenerate biological models)
display('Starting BP continuation');

% set options
opt=contset(opt,'Backward',1);

[x0,v0]=init_BP_BP(@PS_predator_prey,BP,parBP,[ap1 ap2 ap3],2);
[xbp,vbp,sbp,hbp,fbp]=cont(@branchpoint,x0,[],opt); xbp(MM+1,end)
while (xbp(MM+1,end)<5 && xbp(MM+1,end)>1)
    [xbp,vbp,sbp,hbp,fbp]=cont(xbp,vbp,sbp,hbp,fbp,cds); xbp(MM+1,end)
end

% Plot
figure(2)
cpl(xbp,vbp,sbp,[MM+1 MM+2]);

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

% set options
opt=contset(opt,'Backward',0);

[x0,v0]=init_H_H(@PS_predator_prey,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
while (xh(MM+1,end)<5 && xh(MM+1,end)>1)
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
end

% Plot
figure(2)
cpl(xh,vh,sh,[MM+1 MM+2]);

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

% set options
opt=contset(opt,'Backward',1);

[x0,v0]=init_H_H(@PS_predator_prey,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
while (xh(MM+1,end)<5 && xh(MM+1,end)>1)
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
end

% Plot
figure(2)
cpl(xh,vh,sh,[MM+1 MM+2]);
