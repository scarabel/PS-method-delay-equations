% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Diekmann, Vermiglio, Numerical bifurcation analysis of renewal
% equations via pseudospectral approximation, available at 
% https://arxiv.org/abs/2012.05364
  
%% MC_example23.m
% command line instructions for MatCont continuation of the system defined
% in PS_example23.m

clear;
clearvars -global cds

% Discretization parameters
M=20;
    
% Initial parameter values
tau_max=3;
abar=1;
gamma0=-0.2;
par=[gamma0,abar,tau_max,M]';

% Approximated equilibrium corresponding to par
xeq=0;
yeq=[];

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
%ap3=3;
TOL=1e-6;
TestTOL=1e-6;

%% Continuation process

MM=M; % dimension of the approximating ODE system
handles=feval(@PS_specialRE_NBV);
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
opt=contset(opt,'MaxNumPoints',100);

state_eq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_specialRE_NBV,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=0;
while ((length(se)<3) && xe(end,end)< 10 && jj<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
    jj=jj+1;
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

% detection of the equilibrium (or limit cycle) at selected points
UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P ';
opt=contset(opt,'Userfunctions',1);
opt=contset(opt,'UserfunctionsInfo',UserInfo);

% set options
opt=contset(opt,'Backward',0);

[x0,v0]=init_BP_EP(@PS_specialRE_NBV,BP,parBP,sBP,0.001);

tic
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
Time_Eq(M) = toc;

jj=1;
while (xe(end,end)<5 && jj<5)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
    jj=jj+1;
end

% Plot
figure(1)
cpl(xe,ve,se,[MM+1 1]);
xlabel('beta');
title('Bifurcation specialRE')

%%
% Plot of bifurcation diagram of first component
angles = pi*(2*[1:M]'-1)/(2*M); 
Nodes = [0;0.5*tau_max*cos(angles)-0.5*tau_max];
DD = poldif(Nodes,1);
DM = DD(2:end,2:end);

[QuadWeights,QuadNodes]=cheb_quad(50,-tau_max,-abar);

for index_sol = 1:size(xe,2)
    BB=DM*xe(1:end-1,index_sol);
    der = polint(Nodes(2:end),BB,QuadNodes);
    FM = 0.5*exp(xe(end,index_sol))*QuadWeights*(der.*exp(-der));
    b0(index_sol) = FM;
end

figure(10)
plot(xe(end,:),b0,'b');
hold on

for ii=2:length(se)-1
    index=se(ii).index;
    plot(xe(end,index),b0(index),'or');
end

figure
plot(eigs(DM),'o')

% savefig([num2str(M),'_bif']);

% Plot of eigenvalues at selected points
num_hopf=0; % used only to select first Hopf point
for jj=1:size(se,1)
    index = se(jj).index;
    if strcmp(se(jj).label,'P ')==1
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues, M=',num2str(M),', log(gamma)= ',num2str(xe(end,index))]);
%        savefig([num2str(M),'_eig_gamma',num2str(round(xe(end,index)))]);
    elseif (strcmp(se(jj).label,'H ')==1 && num_hopf==0)
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues at H, M=',num2str(M),', gamma= ',num2str(xe(end,index))]);
%        savefig([num2str(M),'_eig_H']);
        num_hopf=1;
    end
end


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

Hopf(M) = xe(end,H_index);

%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

ntst=40; % number of interval
ncol=4; % degree of polynomial

TOL=1e-3;
TestTOL=1e-3;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TestTOL);

% set options
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',1);
opt=contset(opt,'InitStepsize',0.1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'Adapt',0);
opt=contset(opt,'MaxNumPoints',100);

% detection of the limit cycle
UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P ';
opt=contset(opt,'Userfunctions',1);
opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_H_LC(@PS_specialRE_NBV,H,parH,ap1,1e-3,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
jj=0;
while (jj<5 && xlc(end,end)<4.2)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
    jj=jj+1;
end

%% Plot max and min periodic solutions

Per_Solutions = zeros(ntst*ncol+1,size(xlc,2));
for ind_persol=1:size(xlc,2)
    for ind_mesh=1:ntst*ncol+1
        % BB=DM*xlc((ind_mesh-1)*(MM*ncol)+1:(ind_mesh-1)*(MM*ncol)+MM,ind_persol);
        % b0_per(ind_mesh)=polint(tau_max*UnitNodes(2:end),BB,0);
        BB = DM*xlc((ind_mesh-1)*MM+1:(ind_mesh-1)*MM+MM,ind_persol);
        der = polint(Nodes(2:end),BB,QuadNodes);
        b0_per = 0.5*exp(xlc(end,ind_persol))*QuadWeights*(der.*exp(-der));
        Per_Solutions(ind_mesh,ind_persol) = b0_per;

    end

end

upperbound=max(Per_Solutions);
lowerbound=min(Per_Solutions);

figure(10)
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'or',xlc(end,index),lowerbound(index),'or');
end

% savefig([num2str(M),'_bif']);

for jj=1:size(slc,1)
    index = slc(jj).index;
    if strcmp(slc(jj).label,'P ')==1
        figure
        lambda=flc(ntst+2:end,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Multipliers, M=',num2str(M),', log(gamma)= ',num2str(xlc(end,index))]);
%        savefig([num2str(M),'_mult_gamma',num2str(round(xlc(end,index)))]);
    elseif strcmp(slc(jj).label,'PD ')==1
        figure
        lambda=flc(ntst+2:end,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Multipliers at PD, M=',num2str(M),', log(gamma)= ',num2str(xlc(end,index))]);
%        savefig([num2str(M),'_mult_PD']);
    end
end



%% detection of PD bifurcation

SPD=1;

for S=1:size(slc)
    if strcmp(slc(S).label,'PD ')==1
        PD_index=slc(S).index;
        SPD=S;
        break;
    end
end

par_PD=xlc(end,PD_index);
parPD=parH;
parPD(ap1)=par_PD;

Per_Doubl(M) = xlc(end,PD_index);


%% Continuation of periodic solutions from PD
 
% TOL=1e-3;
% TestTOL=1e-3;
% opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
% opt=contset(opt,'TestTolerance',TestTOL);
% 
opt=contset(opt,'InitStepsize',1e-2);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'Multipliers',0);

if M==20
    TestTOL=1e-2;
    opt=contset(opt,'TestTolerance',TestTOL);
end

% opt=contset(opt,'Adapt',0);
% opt=contset(opt,'MaxNumPoints',50);

[xpd0,vpd0]=init_PD_LC(@PS_specialRE_NBV,xlc,slc(SPD),ntst,ncol,1);
[xlc1,vlc1,slc1,hlc1,flc1]= cont(@limitcycle,xpd0,vpd0,opt); xlc1(end,end)
[xlc1,vlc1,slc1,hlc1,flc1]= cont(xlc1,vlc1,slc1,hlc1,flc1,cds); xlc1(end,end)
jj=1;
while (xlc1(end,end)<4.5 && jj<5)
    [xlc1,vlc1,slc1,hlc1,flc1]= cont(xlc1,vlc1,slc1,hlc1,flc1,cds); xlc1(end,end)
    jj=jj+1;
end

%% Plot max and min periodic solutions

ninterp=100;

    % re-interpolation for a smoother plot
    mesh_refined=linspace(0,1,ninterp);
    Per_Solutions = zeros(length(mesh_refined),size(xlc1,2));
    for ind_persol=1:size(xlc1,2)
        Per_Solutions(:,ind_persol) = interp1(flc1(1:ntst+1,ind_persol),xlc1(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
    end
upperbound1=max(Per_Solutions);
lowerbound1=min(Per_Solutions);

figure(1)
plot(xlc1(end,:),upperbound1,'g',xlc1(end,:),lowerbound1,'g');

for ii=2:length(slc1)-1
    index=slc1(ii).index;
    plot(xlc1(end,index),upperbound1(index),'og',xlc1(end,index),lowerbound1(index),'og');
end

% savefig([num2str(M),'_bif']);
% save([num2str(M),'_bif']);



%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

ap2=3;

TOL=1e-6;

% set options
%opt=contset(opt,'MaxStepsize',1e-1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',1);

% % detection of the equilibrium (or limit cycle) at selected points
% opt=contset(opt,'Userfunctions',0);
% opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_H_H(@PS_specialRE_NBV,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
% jj=0;
% while (xh(MM+2,end)>10 && xh(MM+2,end)>0 && jj<5)
%      [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
%       jj=jj+1;
% end

% Plot
figure
cpl(xh,vh,sh,[MM+1 MM+2]); hold on
xlabel('\log\gamma'); ylabel('tau');
title('Hopf curve')


%%
opt=contset(opt,'Backward',1);
[x0,v0]=init_H_H(@PS_specialRE_NBV,H,parH,[ap1 ap2]);
%tic
[xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
%time_hopf=toc

jj=0;
while (xhb(MM+1,end)<10 && jj<10)
    [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
    jj=jj+1;
end

cpl(xhb,vhb,shb,[MM+1 MM+2]); hold on

title(['Regions, M=',num2str(M)]);
% savefig([num2str(M),'_regions']);
% 
% save([num2str(M),'_bif_full']);


% ------------
%%% AUXILIARY FUNCTIONS

function [w,x]=cheb_quad(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
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
