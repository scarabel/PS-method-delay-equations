% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O. (2020)
% Pseudospectral approximation of Hopf bifurcation for delay differential
% equations, SIAM Journal on Applied Dynamical Systems.
%
%% MC_blowflies.m
% MatCont continuation of the Nicholson's blowflies equation
% y'(t) = -mu*y(t)+beta*y(t-1)*exp(-y(t-1)) defined in PS_blowflies

clear;
clearvars -global cds
close all

% Discretization index
M=20;

% Initial parameter values
log_beta=0;
mu=7;
par=[log_beta,mu,M]';

% Approximated equilibrium corresponding to par
yeq=0;

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
TOL=1e-6;
TestTOL=1e-6;

%% Continuation process

MM=M+1; % dimension of the approximating ODE system
handles=feval(@PS_blowflies);
opt=contset;
global cds;

%% Equilibrium continuation from initial point [xeq;yeq]

display('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);

Weq=feval(handles{1},M,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_blowflies,Weq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=1;
while (xe(end,end)<5 &&  jj<5)
     [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
     jj=jj+1;
end

if xe(end,end)<0
    opt=contset(opt,'Backward',1);
    [x0,v0]=init_EP_EP(@PS_blowflies,Weq,par0,ap1);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    jj=1;
    while (xe(end,end)<5 &&  jj<5)
        [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
        jj=jj+1;
    end
    
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
hold on;
xlabel('log(beta)')
%axis([0 10 0 1.4])

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
opt=contset(opt,'Backward',0);
%opt=contset(opt,'MaxStepsize',1e-1);
%opt=contset(opt,'MaxNumPoints',10);

% detection of the equilibrium point
UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P ';
opt=contset(opt,'Userfunctions',1);
opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_BP_EP(@PS_blowflies,BP,parBP,sBP,0.01);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=1;
while (xe(end,end)<5 &&  jj<5)
     [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
     jj=jj+1;
end

if xe(end,end)<0
    opt=contset(opt,'Backward',1);
    [x0,v0]=init_EP_EP(@PS_blowflies,Weq,par0,ap1);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    jj=1;
    while (xe(end,end)<5 &&  jj<5)
        [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
        jj=jj+1;
    end

end

figure(1);
cpl(xe,ve,se,[MM+1 1]);

%% Plot of eigenvalues in three different points
num_hopf=0;
for jj=1:size(se,1)
    index = se(jj).index;
    if strcmp(se(jj).label,'P ')==1
        ratio = exp(xe(end,index))/mu;
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues, M=',num2str(M),', beta/mu= ',num2str(ratio), ', mu=',num2str(mu)]);
%         savefig([num2str(M),'_eig',num2str(jj), ', mu=',num2str(mu)]);
                
        % plot points in two-parameter plane
        figure(9)
        plot(mu,ratio,'x'); hold on
    elseif (strcmp(se(jj).label,'H ')==1 && num_hopf==0)
        ratio = exp(xe(end,index))/mu;
        num_hopf=1;
        figure
        lambda=fe(:,index);
        plot(real(lambda),imag(lambda),'o');
        title(['Eigenvalues at Hopf, M=',num2str(M),', beta/mu= ',num2str(ratio), ', mu=',num2str(mu)]);
%         savefig([num2str(M),'_eig_Hopf', ', mu=',num2str(mu)]);
        
        % plot points in two-parameter plane
        figure(9)
        plot(mu,ratio,'x'); hold on
    end
end

figure(9)
xlabel('mu'); ylabel('beta/mu'); 
title('Points for computation of eigenvalues')
% savefig([num2str(M),'_Points_eigs'])

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

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

ap2=2;

%TOL=1e-2;
% set options
% opt=contset(opt,'MaxStepsize',1e-1);
% opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
% opt=contset(opt,'TestTolerance',TOL);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxNumPoints',200);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'Userfunctions',0);

[x0,v0]=init_H_H(@PS_blowflies,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+2,end)
jj=0;
while (xh(MM+1,end)<10 && jj<5)
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
      jj=jj+1;
end
 
% Plot
figure(9); hold on
plot(xh(MM+2,:), exp(xh(MM+1,:))./xh(MM+2,:))
axis([0 5 0 50])
xlabel('mu'); ylabel('beta/mu');
title(['Hopf curve approximated with Matcont, M=',num2str(M)])

% savefig([num2str(M),'_Hopf_curve'])

% Plot imaginary part of eigenvalue
figure(10); clf
Omega=imag(fh(end-1,:));
plot(xh(MM+2,:), Omega);
xlabel('mu'); ylabel('omega'); title(['Imaginary part of eigenvalue, M=',num2str(M)])
% figname=[num2str(M),'_imag_part',];
% savefig(figname);

% Plot Lyapunov coefficient
figure(11)
plot(xh(MM+2,:),hh(6,:))
xlabel('mu'); ylabel('l1'); title(['First Lyapunov coefficient, M=',num2str(M)])
% figname=[num2str(M),'_1lyapunov'];
% savefig(figname);

% omega = fh(end-1,:);
% 
% % normalization
% [pesi,Nodes,D,q]=cheb_delay(M,-1,0);
% normalization = abs(sum(exp(2*Nodes.*omega),1));
% for jj= 1: length(omega)
%     normalization(jj)= exp(omega(jj)*Nodes')*exp(omega(jj)*Nodes')';
% end
% 
% l1 = hh(6,:).*normalization;
% %l1_bis = hh(6,:).*normalization./abs(omega);
% 
% figure
% plot(xh(MM+2,:),l1)
% xlabel('mu'); ylabel('l1'); title(['First Lyapunov coefficient (scaled as DDE), M=',num2str(M)])
% figname=['1lyapunov_scaled_',num2str(M)];
% savefig(figname);

%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

ntst=40; % number of interval
ncol=4; % degree of polynomial

%TOL=1e-3;
%TestTOL=1e-3;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TestTOL);

% set options
opt=contset(opt,'Backward',0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'MaxStepsize',1);
%opt=contset(opt,'MinStepsize',1e-3);
%opt=contset(opt,'InitStepsize',1e-3);
%opt=contset(opt,'Adapt',0);

% detection of the equilibrium point
UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P ';
opt=contset(opt,'Userfunctions',1);
opt=contset(opt,'UserfunctionsInfo',UserInfo);

[x0,v0]=init_H_LC(@PS_blowflies,H,parH,ap1,1e-3,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
jj=0;
while (xlc(end,end)<5 && jj<10)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
    jj=jj+1;
end

%% Plot max and min periodic solutions

ninterp=100;

    % re-interpolation for a smoother plot
    mesh_refined=linspace(0,1,ninterp);
    Per_Solutions = zeros(length(mesh_refined),size(xlc,2));
    for ind_persol=1:size(xlc,2)
        Per_Solutions(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
    end
upperbound=max(Per_Solutions);
lowerbound=min(Per_Solutions);

figure(1); hold on
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'og',xlc(end,index),lowerbound(index),'og');
end


% re-interpolation for a smoother plot
ninterp=100;
mesh_refined=linspace(0,1,ninterp);

Per_Solution = zeros(length(mesh_refined),1);

for jj=1:size(slc,1)
    index = slc(jj).index;
    if strcmp(slc(jj).label,'P ')==1
        ind_persol= slc(jj).index;
        T = xlc(end-1,index); 
        Per_Solution(:) = interp1(flc(1:ntst+1,ind_persol),xlc(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        figure(10);
        plot([T*mesh_refined,T*(1+mesh_refined)],[Per_Solution;Per_Solution]);
        title(['Two periods, log beta=',num2str(xlc(end,ind_persol)), ', period=',num2str(xlc(end-1,ind_persol))])
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

% %% Continuation of periodic solutions from PD
%  
% TOL=1e-3;
% TestTOL=1e-3;
% opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
% opt=contset(opt,'TestTolerance',TestTOL);
% 
% % opt=contset(opt,'InitStepsize',1e-4);
% %opt=contset(opt,'MaxStepsize',1e-1);
% %opt=contset(opt,'Adapt',0);
% %opt=contset(opt,'MaxNumPoints',20);
% 
% [xpd0,vpd0]=init_PD_LC(@PS_blowflies,xlc,slc(SPD),ntst,ncol,1);
% [xlc1,vlc1,slc1,hlc1,flc1]= cont(@limitcycle,xpd0,vpd0,opt); xlc1(end,end)
% jj=1;
% while (xlc1(end,end)<6 && jj<5)
%     [xlc1,vlc1,slc1,hlc1,flc1]= cont(xlc1,vlc1,slc1,hlc1,flc1,cds); xlc1(end,end)
%     jj=jj+1;
% end
% 
% %% Plot max and min periodic solutions
% 
% ninterp=100;
% 
%     % re-interpolation for a smoother plot
%     mesh_refined=linspace(0,1,ninterp);
%     Per_Solutions = zeros(length(mesh_refined),size(xlc1,2));
%     for ind_persol=1:size(xlc1,2)
%         Per_Solutions(:,ind_persol) = interp1(flc1(1:ntst+1,ind_persol),xlc1(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
%     end
% upperbound1=max(Per_Solutions);
% lowerbound1=min(Per_Solutions);
% 
% %figure; hold on
% plot(xlc1(end,:),upperbound1,'r',xlc1(end,:),lowerbound1,'r');
% 
% for ii=2:length(slc1)-1
%     index=slc1(ii).index;
%     plot(xlc1(end,index),upperbound1(index),'og',xlc1(end,index),lowerbound1(index),'og');
% end
% 
% 
% for jj=1:size(slc1,1)
%     index = slc1(jj).index;
%     if strcmp(slc1(jj).label,'P ')==1
%         ind_persol= slc1(jj).index;
%         T = xlc1(end-1,index); 
%         Per_Solution(:) = interp1(flc1(1:ntst+1,ind_persol),xlc1(1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
%         figure(10); subplot(2,1,2)
%         plot([T*mesh_refined,T*(1+mesh_refined)],[Per_Solution;Per_Solution]);
%         title(['Two periods, log beta=',num2str(xlc1(end,ind_persol)), ', period=',num2str(xlc1(end-1,ind_persol))])
%         savefig([num2str(M),'_per_sol_ntst',num2str(ntst)])
%     end
% end


%% PD continuation in two parameters % ap1,ap2 = index of continuation
% parameters in the vector par display('Starting PD continuation');

ap2=2;

TOL=1e-3; % set options opt=contset(opt,'MaxStepsize',1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL); 
opt=contset(opt,'Singularities',0);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Multipliers',0); 
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'MaxStepsize',0.5);
%opt=contset(opt,'InitStepsize',0.01);

[x0,v0]=init_PD_PD(@PS_blowflies,xlc,slc(SPD),[ap2 ap1],ntst,ncol);
[xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); % 
[xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end,end) %       
[xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end,end) %       
l=size(xpd,1);
% jj=1;
% while (xpd(l-1,end)<5 && jj<10)
%     [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(l-1,end) %       
%     jj=jj+1;
% end

%%
% opt=contset(opt,'Backward',1);
% [x0,v0]=init_PD_PD(@PS_blowflies,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
% [xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); % 
% [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end,end) %       

% jj=0; 
% while (xpd(end,end)>2.65 && jj<5) %
% [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end,end) %       
% jj=jj+1;
% end %

% %%
% opt=contset(opt,'Backward',1);
% [x0,v0]=init_PD_PD(@PS_blowflies,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
% [xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); % 
% jj=0; 
% while (xpd(end,end)>2.65 && jj<5) %
% [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end,end) %       
% jj=jj+1;
% end %

%% Plot of stability regions


figure
l=size(xpd,1);
%cpl(xpd,vpd,spd,[l l-1]); hold on
plot(xpd(l-1,1:20:end),exp(xpd(l,1:20:end))./xpd(l-1,1:20:end),'r');
xlabel('mu'); ylabel('beta/mu');
title('PD curve in two parameter plane')
% savefig([num2str(M),'_PD_curve_ntst_',num2str(ntst)]);

% % Plot % figure(2) % cpl(xh,vh,sh,[MM+2 MM+1]); hold on % xlabel('mu'); ylabel('beta');


%% Plot of stability regions
figure
plot([0;5], [1;1], 'g'); hold on % existence
plot(xh(MM+2,:),xh(MM+1,:)./xh(MM+2,:),'b'); % Hopf
plot(xpd(end,:),xpd(end-1,:)./xpd(end,:),'r'); % Period doubling
axis([0 5 0 30]);
xlabel('mu'); ylabel('beta/mu');
title(['Stability regions, M=',num2str(M)]);

legend('Existence','Hopf','PD');

% figname=[num2str(M),'_regions'];
% savefig(figname);

% namefile = [num2str(M),'_DDE_blowflies_ntst_',num2str(ntst)];
% save(namefile)
