% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
%
%% Building blocks Matcont continuation
% systemfile.m = Matcont system definition file
% d1,d2 = dimension of RE and DDE, respectively
% M = discretization index
% MM = d1*M+d2*(M+1) = total system dimension

handles=feval(@systemfile);
opt=contset; 
global cds

%% Equilibrium continuation from initial point [xeq;yeq]
% xeq = equilibrium of RE
% yeq = equilibrium of DDE
% par = vectors of parameters
% ap1 = iindex of the continuation parameter in par

disp('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',500);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',1);

Weq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@systemfile,Weq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< UpperBound)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
end

% Plot
% figure(1); clf;
% cpl(xe,ve,se,[MM+1 1]);
% hold on
% or equivalently
% plot(xe(:,end),xe(:,1))


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

% LP, limit point

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

%% Equilibrium continuation from BP
% xeBP = vector of variables at BP
% parBP = parameter vector at BP
% sBP = information about BP from previous continuation
disp('Starting equilibrium continuation from BP');

% set options
% opt=contset(opt,'MaxNumPoints',50);

[x0,v0]=init_BP_EP(@systemfile,BP,parBP,sBP,0.1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< UpperBound)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
end

% Plot
% figure(1)
% cpl(xe,ve,se,[MM+1 1]);
% or equivalently
% plot(xe(:,end),xe(:,1))


%% Limit cycle continuation from H
% xeH = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
disp('Starting LC continuation from H');

% set options
% opt=contset(opt,'MaxNumPoints',50);

ntst=20; % number of intervals
ncol=4; % degree of polynomials

[x0,v0]=init_H_LC(@systemfile,H,parH,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< UpperBound)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
% upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
% lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
% figure(1)
% plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
% xlabel('parameter','interpreter','latex');
% ylabel('max/min','interpreter','latex')
% or equivalently my own function cpl_per:
% cpl_per(xlc,slc,flc,ntst,ncol,MM,M+1,size(xlc,2));

%% Plot profile of orbit corresponding to index ind_persol
ind_persol=size(xlc,2);
%ind_persol=1;
p_persol=xlc(end,ind_persol);
period_persol=xlc(end-1,ind_persol);
mesh_persol= xlc(end-1,end)*flc(1:ntst+1,ind_persol);
profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol);

% Compute the delay along the periodic solution:
delay_persol=zeros(ncol*ntst+1,1);
for index_mesh=1:ntst*ncol+1
    state_mesh=xlc((index_mesh-1)*MM+1:(index_mesh*MM),end);
    delay_persol(index_mesh)= feval(handles{10},state_mesh,a,p_persol,ka,kp,kg,mu,muw,eps,alpha,x1,x2,aux,M);
end

tt_persol=linspace(0,period_persol,1000);
sol=interp1(mesh_persol,profile,tt_persol);

figure
plot(mesh_persol,profile,'.'); hold on
plot(tt_persol,sol,'b-');
plot(mesh_persol,eq_ref(1),'g')
xlabel('time','interpreter','latex');
ylabel('$state$','interpreter','latex');

%% Limit cycle continuation from PD
% xlcPD = output of the continuation where PD was detected
% parPD = parameter vector at PD
% ap1 = index of continuation parameter in vector par
disp('Starting LC continuation from PD');

% set options
% opt=contset(opt,'MaxNumPoints',50);

ntst=20; % numer of interval
ncol=4; % degree of polynomial

[x0,v0]=init_PD_LC(@systemfile,xlcPD,sPD,ntst,ncol,1e-6);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
while ((length(slc)<3) && xlc(end,end)< UpperBound)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
end

% Plot max and min periodic solutions
% upperbound=max(xlc(1:MM:((ntst*ncol+1)*MM),:));
% lowerbound=min(xlc(1:MM:((ntst*ncol+1)*MM),:));
% figure(1)
% plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
% xlabel('parameter','interpreter','latex');
% ylabel('max/min','interpreter','latex')
% or equivalently my own function cpl_per:
% cpl_per(xlc,slc,flc,ntst,ncol,MM,M+1,size(xlc,2));

% Plot profile of orbit corresponding to index j
% j='end';
% figure
% mesh= xlc(end-1,j)*flc(1:ntst+1,j);
% profile= xlc(1:MM*ncol:((ntst*ncol+1)*MM),j);
% plot(mesh,profile);

%% BP continuation in two parameters
% BP = vector of variables at BP
% parBP = parameter vector at BP
% ap1,ap2,ap3 = index of continuation parameter in par (often ap3 is
% auxiliary parameter in degenerate biological models)
disp('Starting BP continuation');

% set options
% opt=contset(opt,'MaxNumPoints',50);

[x0,v0]=init_BP_BP(@systemfile,BP,parBP,[ap1 ap2 ap3],1);
[xbp,vbp,sbp,hbp,fbp]=cont(@branchpoint,x0,v0,opt);
for j=1:2
    [xbp,vbp,sbp,hbp,fbp]=cont(xbp,vbp,sbp,hbp,fbp,cds); xbp(MM+1,end)
end

% Plot
% figure(2); clf;
% cpl(xbp,vbp,sbp,[MM+1 MM+2]);
% hold on
% or equivalently:
% plot(xbp(MM+1,:),xbp(MM+2,:))


%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
disp('Starting H continuation');

% set options
% opt=contset(opt,'MaxNumPoints',50);

[x0,v0]=init_H_H(@systemfile,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt);
for j=1:3
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
end

% Plot
% figure(2)
% cpl(xh,vh,sh,[MM+1 MM+2]);
% or equivalently:
% plot(xh(MM+1,:),xh(MM+2,:))
