% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Diekmann, Vermiglio, Numerical bifurcation analysis of renewal
% equations via pseudospectral approximation, available at 
% https://arxiv.org/abs/2012.05364

% MC_example21
% command line instructions for MatCont continuation of the system defined
% in PS_example21.m

clear;
clearvars -global cds
close all

% Discretization parameters
M=40;

% Initial parameter values
d1=1; d2=0; % dimension of RE and DDE, respectively
loggamma=-1;
q=0.6;
tau=1;

par=[loggamma,q,tau,M]';

%% Computation of Fourier coefficients

k=3; theta=0.1; 
kernel_shape = @(s) s.^(k-1).*exp(-s./theta).*(s<=tau);
kernel = @(s) kernel_shape(s)/integral(@(sigma) kernel_shape(sigma),0,tau);
figure(23)
plot(0:0.1:tau, kernel(0:0.1:tau),'r')

n_quad = M-1; %2*M; %max([M,20]);
[QuadWeights,QuadNodes,~,~]=cheb_delay(n_quad,-tau,0);
kernel_shape = @(s) s.^(k-1).*exp(-s./theta).*(s<=tau);
norm_const = QuadWeights*kernel_shape(-QuadNodes); %integral(@(sigma) kernel_shape(sigma),0,tau);
kernel = @(s) kernel_shape(s)/norm_const;
hold on
plot(0:0.1:tau, kernel(0:0.1:tau),'b')

%figure(10); clf
%plot(0:0.01:1, kernel(0:0.01:1));

for nn=1:5
    Fourier(nn) = 2*integral(@(s) kernel(s).*sin(2*pi*nn*s),0,1)
end

%%
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
handles=feval(@PS_prelude);
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
opt=contset(opt,'MaxNumPoints',100);
% opt=contset(opt,'MaxStepsize',1e-1);

opt=contset(opt,'Backward',0);
state_eq=feval(handles{1},M,xeq,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_prelude,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
while ((length(se)<3) && xe(end,end)< 10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
hold on;

%% Detection of singular points
% xe,ve,se,he,fe = output of previous continuation
% par = current parameter vector
% ap1 = index of continuation parameter in vector par

% BP, branching point
BP_index=[];
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
opt=contset(opt,'MaxNumPoints',500);

opt=contset(opt,'Backward',0);
% 
% % detection of the equilibrium (or limit cycle) at selected points
% UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P ';
% opt=contset(opt,'Userfunctions',1);
% opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_BP_EP(@PS_prelude,BP,parBP,sBP,0.001);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=1;
while ((length(se)<3) && xe(end,end)<4 && jj<5)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
    jj=jj+1;
end

if xe(end,end)<0
    opt=contset(opt,'Backward',1);
    [x0,v0]=init_BP_EP(@PS_prelude,BP,parBP,sBP,0.001);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    jj=1;
    while ((length(se)<3) && xe(end,end)<10 && jj<5)
        [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
        jj=jj+1;
    end
    opt=contset(opt,'Backward',0);
end

%% Plot
figure(1)
cpl(xe,ve,se,[MM+1 1]);
xlabel('beta');
title('Bifurcation prelude')

% Plot of bifurcation diagram of first component
angles = pi*(2*[1:M]'-1)/(2*M); 
Nodes = [0;0.5*tau*cos(angles)-0.5*tau];
DD = poldif(Nodes,1);
DM = DD(2:end,2:end);

p = pi*(2*(0:M-1)'+1)/(2*M);
x=[1;sin(pi/2-p)]; % nodes with addition of 1 % either cos(p) or sin (pi/2-p) 
X=repmat(x,1,M+1);
dX=X-X';
c=[2^(M-1)/M*prod(dX(1,2:end)); ((-1).^(0:M-1)').*dX(2:end,1)./sin(p)];
D=(c*(1./c'))./(dX+(eye(M+1)));
D=D-diag(sum(D')); % differentiation matrix
% scaling
Nodes = 0.5*tau*(x-1);
DM = 2/tau*D(2:end,2:end);

n_quad = M-1;
[QuadWeights,QuadNodes]=cheb_quad(nquad,-tau,0);

b0 = zeros(size(xe,2),1);
for index_sol = 1:size(xe,2)
   
    derb = DM*xe(1:end-1,index_sol);
    der_quad = polint(Nodes(2:end),derb,QuadNodes);
    int1 = QuadWeights*(kernel(-QuadNodes).*der_quad);

    FM = exp(xe(end,index_sol))*(1-QuadWeights*der_quad)*int1;
    b0(index_sol) = FM; % select component in node -tau
end

figure(10)
plot(xe(end,:),b0,'b');
hold on

for ii=2:length(se)-1
    index=se(ii).index;
    plot(xe(end,index),b0(index),'or');
end

% Plot of eigenvalues at selected points
num_hopf=0; % used only to select first Hopf point
for jj=1:size(se,1)
    index = se(jj).index;
    if strcmp(se(jj).label,'P ')==1
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues, M=',num2str(M),', log(gamma)= ',num2str(xe(end,index))]);
        % savefig([num2str(M),'_eig_gamma',num2str(round(xe(end,index)))]);
    elseif (strcmp(se(jj).label,'H ')==1)
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues at H, M=',num2str(M),', gamma= ',num2str(xe(end,index))]);
        % savefig([num2str(M),'_eig_H_',num2str(num_hopf)]);
        num_hopf=num_hopf+1;
    end
end

%% Detection of singular points
% H, Hopf point

num_hopf=0;
H_index=[];
for ii=1:size(se)
    if strcmp(se(ii).label,'H ')==1
        num_hopf=num_hopf+1;
%        Hopf(M,num_hopf) = xe(end,se(ii).index);
        if num_hopf==1
            H_index=se(ii).index;
        end
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
TestTOL=1e-6; TOL=1e-6;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TestTOL);

opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepSize',0.1);
%opt=contset(opt,'MaxStepSize',5e-2); % M=15
%opt=contset(opt,'Adapt',3);

ntst=40; % number of intervals
ncol=4; % degree of polynomial

% % detection of the equilibrium (or limit cycle) at selected points
% opt=contset(opt,'Userfunctions',1);
% opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_H_LC(@PS_prelude,H,parH,ap1,0.1,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
% jj=0;
% while (xlc(end,end)<4 && jj<10)
%     [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
%     jj=jj+1;
% end

% save([num2str(M),'_prelude_ntst_',num2str(ntst)]);

%% Plot max and min periodic solutions
Per_Solutions = zeros(ntst*ncol+1,size(xlc,2));
for ind_persol=1:size(xlc,2)
    for ind_mesh=1:ntst*ncol+1
        % BB=DM*xlc((ind_mesh-1)*(MM*ncol)+1:(ind_mesh-1)*(MM*ncol)+MM,ind_persol);
        % b0_per(ind_mesh)=polint(tau_max*UnitNodes(2:end),BB,0);
        derb = DM*xlc((ind_mesh-1)*MM+1:(ind_mesh-1)*MM+MM,ind_persol);
        der_quad = polint(Nodes(2:end),derb,QuadNodes);
        int1 = QuadWeights*(kernel(-QuadNodes).*der_quad);
        FM = exp(xlc(end,ind_persol))*(1-QuadWeights*der_quad)*int1;
        b0_per = FM;
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

% savefig([num2str(M),'_bif_ntst_',num2str(ntst)]);

%% save one periodic solution
index_per = 83;
for jj=1:ntst
Mesh(ncol*(jj-1)+1:ncol*jj+1) = linspace(flc(jj,index_per),flc(jj+1,index_per),ncol+1);
end
figure
plot(Mesh,Per_Solutions(:,index_per))

%%

ind_PD=1;
for jj=1:size(slc,1)
    index = slc(jj).index;
    if strcmp(slc(jj).label,'P ')==1
        figure
        lambda=flc(ntst+2:end,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Multipliers, M=',num2str(M),', log(gamma)= ',num2str(xlc(end,index))]);
        % savefig([num2str(M),'_mult_gamma',num2str(round(xlc(end,index)))]);
    elseif strcmp(slc(jj).label,'PD ')==1 && ind_PD==1
        figure
        lambda=flc(ntst+2:end,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Multipliers at PD, M=',num2str(M),', gamma= ',num2str(xlc(end,index))]);
        % savefig([num2str(M),'_mult_PD_ntst_',num2str(ntst)]);
        ind_PD=ind_PD+1;
    end
end



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
