% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O. (2020)
% Pseudospectral approximation of Hopf bifurcation for delay differential
% equations, SIAM Journal on Applied Dynamical Systems.
%
%% systems_ddebiftool.m
% Bifurcation analysis of the neural system using dde-biftool
% The code requires that the folder containing dde-biftool codes is on the
% Matlab path. 

clear
close all

% DEFINITIONS
% model definition
% p = [k,c,tau]
k0=2; 
c0=1;
tau0 = 1;
par0 = [k0, c0, tau0];

x0 = [c0;2/(k0*c0^2)];

ind_k = 1;
ind_c = 2;
ind_tau = 3;

system_def = @(x,p) [1-0.5*p(1)*x(1,1)*x(1,2)*x(2,2); x(1,1)-p(2)];

funcs = set_funcs(...
    'sys_rhs',system_def,...
    'sys_tau', @() ind_tau);

% choice of continuation parameter
ind_contpar = ind_k;

%% Initialization of branch of non-trivial equilibria

maxstep=0.1; maxbound=6; 

nontriv_eqs=SetupStst(funcs,'x',x0,'parameter',par0,'step',0.1,...
    'contpar',ind_contpar,'max_step',[ind_contpar,maxstep],'max_bound',[ind_contpar,maxbound]);

%% Compute and find stability of non-trivial equilibria

figure(1);clf
nontriv_eqs=br_contn(funcs,nontriv_eqs,100);
%nontriv_eqs=br_rvers(nontriv_eqs);
%nontriv_eqs=br_contn(funcs,nontriv_eqs,100);

nontriv_eqs=br_stabl(funcs,nontriv_eqs,0,1);
nunst_eqs=GetStability(nontriv_eqs);
ind_hopf=find(nunst_eqs<2,1,'last');
fprintf('Hopf bifurcation near point %d\n',ind_hopf);

%% Continue Hopf bifurcation in two parameters
[hbranch,suc]=SetupHopf(funcs,nontriv_eqs,ind_hopf,...
    'contpar',[ind_k,ind_c],'dir',ind_k,'step',0.1,'max_bound',[ind_c,maxbound]);
hbranch.method.continuation.plot=1;

figure(2);clf
hbranch=br_contn(funcs,hbranch,100);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,300);
xlabel('k'); ylabel('c');

for jj = 1:size(hbranch.point,2)
    k_hopf(jj) = hbranch.point(jj).parameter(ind_k);
    c_hopf(jj) = hbranch.point(jj).parameter(ind_c);
end

figure(5); clf
plot(k_hopf, c_hopf); hold on
axis([0 5 0 3])
xlabel('k'); ylabel('c')
% savefig('db_hopf');



%% Branch off at  Hopf bifurcation
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');

maxbound=10;

[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',1,'intervals',40,'degree',4,...
    'max_bound',[ind_contpar,maxbound],'max_step',[ind_contpar,maxstep],'radius',0.1,'step',0.1);
if ~suc
    error('fail',...
        ' initialization of periodic orbit failed');
end

figure(1);
hold on
per_orb=br_contn(funcs,per_orb,100);
per_orb=br_stabl(funcs,per_orb,0,1);
nunst_per=GetStability(per_orb,'exclude_trivial',true);


%% Find period doubling bifurcations in two parameters
ind_pd=find(diff(nunst_per)==1);
[pdfuncs,pdbranch1,suc]=SetupPeriodDoubling(funcs,per_orb,ind_pd(1),...
    'contpar',[ind_k,ind_c],'dir',ind_k,'step',1e-2);
if ~suc
    error('fail',...
        ' initialization of period doubling failed');
end
figure(2);
pdbranch1=br_contn(pdfuncs,pdbranch1,100);
pdbranch=br_rvers(pdbranch1);
pdbranch=br_contn(pdfuncs,pdbranch,300);


for jj = 1:size(pdbranch.point,2)
    k_pd(jj) = pdbranch.point(jj).parameter(ind_k);
    c_pd(jj) = pdbranch.point(jj).parameter(ind_c);
end

figure(5); hold on
plot(k_pd, c_pd); hold on
axis([0 5 0 3])
xlabel('k'); ylabel('c')
% savefig('db_regions');


