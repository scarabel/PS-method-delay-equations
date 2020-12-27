% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O. (2020)
% Pseudospectral approximation of Hopf bifurcation for delay differential
% equations, SIAM Journal on Applied Dynamical Systems.
%
%% MC_system.m
% MatCont continuation of the neural system
% w'(t) = 1-0.5*k*w(t)*w(t-1)*q(t-1)
% q'(t) = w(t)-c
% Using Chebyshev zeros (plus 0 and -1)
% The code uses the code poldif.m from the Differentiation Matrix Suite
% (Weideman, Reddy, 2000)

clear;
clearvars -global cds
close all

% Discretization index
M=20;

% Initial parameter values
k=0.1; 
c=3;
par=[k,c,M]';

% Approximated equilibrium corresponding to par
yeq=[c;2/(k*c^2)];

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
%ap2=2;
TOL=1e-6;
TestTOL=1e-6;

%% Continuation process

MM=2*(M+1); % dimension of the approximating ODE system
handles=feval(@PS_system);
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
opt=contset(opt,'Backward',1);

Weq=feval(handles{1},M,yeq); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_system,Weq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=1;
while (xe(end,end)<5 &&  jj<5)
     [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
     jj=jj+1;
end

figure(1); clf;
cpl(xe,ve,se,[MM+1 1]);
hold on;
xlabel('k')
%axis([0 10 0 1.4])

% if xe(end,end)<0
%     opt=contset(opt,'Backward',1);
%     [x0,v0]=init_EP_EP(@PS_system,Weq,par0,ap1);
%     [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
%     jj=1;
%     while (xe(end,end)<5 &&  jj<5)
%         [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
%         jj=jj+1;
%     end
%     figure(1); clf;
%     cpl(xe,ve,se,[MM+1 1]);
%     hold on;
%     xlabel('log(beta)')
%     %axis([0 10 0 1.4])
% end


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
opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxNumPoints',200);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);

[x0,v0]=init_H_H(@PS_system,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
% [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)

opt=contset(opt,'Backward',1);
[x0,v0]=init_H_H(@PS_system,H,parH,[ap1 ap2]);
[xh1,vh1,sh1,hh1,fh1]=cont(@hopf,x0,[],opt); xh1(MM+1,end)
%[xh1,vh1,sh1,hh1,fh1]=cont(xh1,vh1,sh1,hh1,fh1,cds); xh1(MM+1,end)
 
% Plot
figure(2); hold on
plot(xh(MM+1,:), xh(MM+2,:))
plot(xh1(MM+1,:), xh1(MM+2,:))
%plot(xh(MM+2,:), exp(xh(MM+1,:))./xh(MM+2,:))
axis([0 5 0 3])
xlabel('k'); ylabel('c');
title(['Hopf curve approximated with Matcont, M=',num2str(M)])

%savefig([num2str(M),'_Hopf_curve'])

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

[x0,v0]=init_H_LC(@PS_system,H,parH,ap1,1e-3,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
jj=0;
while (xlc(end,end)<5 && jj<5)
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
%xlabel('k','interpreter','latex');
%ylabel('max/min','interpreter','latex')

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'og',xlc(end,index),lowerbound(index),'og');
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
% [xpd0,vpd0]=init_PD_LC(@PS_system,xlc,slc(SPD),ntst,ncol,1);
% [xlc1,vlc1,slc1,hlc1,flc1]= cont(@limitcycle,xpd0,vpd0,opt); xlc1(end,end)
% [xlc1,vlc1,slc1,hlc1,flc1]= cont(xlc1,vlc1,slc1,hlc1,flc1,cds); xlc1(end,end)
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
%     end
% end
% 
% 
% %% PD continuation in two parameters % ap1,ap2 = index of continuation
% % parameters in the vector par display('Starting PD continuation');
% 
% ap2=2;
% 
% TOL=1e-3; % set options opt=contset(opt,'MaxStepsize',1);
% opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
% opt=contset(opt,'TestTolerance',TOL); 
% opt=contset(opt,'Singularities',0);
% opt=contset(opt,'Eigenvalues',0);
% opt=contset(opt,'Multipliers',0); 
% opt=contset(opt,'Backward',1);
% opt=contset(opt,'MaxNumPoints',50);
% opt=contset(opt,'MaxStepsize',1);
% 
% [x0,v0]=init_PD_PD(@PS_system,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
% [xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); xpd(end-1,end)  % 
% jj=1;
% while (jj<5 && xpd(end-1,end)<5 && xpd(end,end)<5)
%     [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); xpd(end-1,end) %       
%     jj=jj+1;
% end

% %%
% opt=contset(opt,'Backward',1);
% [x0,v0]=init_PD_PD(@PS_system,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
% [xpd1,vpd1,spd1,hpd1,fpd1]=cont(@perioddoubling,x0,v0,opt); xpd1(end-1,end) %    
% jj=1;
% while (jj<5 && xpd1(end-1,end)<5 && xpd1(end,end)<5)
%     [xpd1,vpd1,spd1,hpd1,fpd1]=cont(xpd1,vpd1,spd1,hpd1,fpd1,cds); xpd1(end-1,end) %   
%     jj=jj+1;
% end

%%

figure
l=size(xpd,1);
%cpl(xpd,vpd,spd,[l l-1]); hold on
%plot(xpd(l-1,:),exp(xpd(l,:))./xpd(l-1,:),'r');
plot(xpd(l-1,:),xpd(l,:),'r'); hold on
%plot(xpd1(l-1,:),xpd1(l,:),'r');
xlabel('k'); ylabel('c');
title('PD curve in two parameter plane')
% savefig([num2str(M),'_PD_curve_ntst',num2str(ntst)]);

% namefile = [num2str(M),'_DDE_system'];
% save(namefile)

%% Plot of stability regions
figure
plot(xh(MM+1,:),xh(MM+2,:),'b'); hold on % Hopf
plot(xh1(MM+1,:),xh1(MM+2,:),'b'); % Hopf
plot(xpd(l-1,:),xpd(l,:),'r.'); % Period doubling
%plot(xpd1(l-1,:),xpd1(l,:),'r'); % Period doubling
axis([0 5 0 5]);
xlabel('k'); ylabel('c');
title(['Stability regions, M=',num2str(M),'_ntst',num2str(ntst)]);

legend('Hopf','PD');

% figname=[num2str(M),'_regions_ntst',num2str(ntst)];
% savefig(figname);
% 
% namefile = [num2str(M),'_DDE_system'];
% save(namefile)

%% Plot periodic solutions

% re-interpolation for a smoother plot
ninterp=100;
mesh_refined=linspace(0,1,ninterp);

Per_Solution = zeros(length(mesh_refined),1);

% first component
ind_comp=1;
for jj=1:size(slc1,1)
    index = slc1(jj).index;
    if strcmp(slc1(jj).label,'P ')==1
        ind_persol= slc1(jj).index;
        T = xlc1(end-1,index); 
        Per_Solution(:) = interp1(flc1(1:ntst+1,ind_persol),xlc1(ind_comp:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        figure(11); subplot(2,1,1)
        plot([T*mesh_refined,T*(1+mesh_refined)],[Per_Solution;Per_Solution]);
        title(['Two periods, k=',num2str(xlc1(end,ind_persol)), ', period=',num2str(xlc1(end-1,ind_persol)),', ind comp=',num2str(ind_comp)])
%         savefig([num2str(M),'_per_sol_ntst',num2str(ntst),'_st'])
        break;
    end
end

% second component
ind_comp=2;
for jj=1:size(slc1,1)
    index = slc1(jj).index;
    if strcmp(slc1(jj).label,'P ')==1
        ind_persol= slc1(jj).index;
        T = xlc1(end-1,index); 
        Per_Solution(:) = interp1(flc1(1:ntst+1,ind_persol),xlc1(ind_comp:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        figure(11); subplot(2,1,2)
        plot([T*mesh_refined,T*(1+mesh_refined)],[Per_Solution;Per_Solution]);
        title(['Two periods, k=',num2str(xlc1(end,ind_persol)), ', period=',num2str(xlc1(end-1,ind_persol)),', ind comp=',num2str(ind_comp)])
        break;
    end
end
