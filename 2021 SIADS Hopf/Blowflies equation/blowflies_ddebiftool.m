% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% De Wolff B, Scarabel F, Verduyn Lunel S, Diekmann O. (2020)
% Pseudospectral approximation of Hopf bifurcation for delay differential
% equations, SIAM Journal on Applied Dynamical Systems.
%
%% blowflies_ddebiftool.m
% Bifurcation analysis of Nicholson's blowflies equation using dde-biftool
% The code requires that the folder containing dde-biftool codes is on the
% Matlab path. 

clear
close all

% DEFINITIONS
% model definition
% p = [log_beta, mu]
log_beta0 = 3;
mu0 = 7;
tau0 = 1;
par0 = [log_beta0, mu0, tau0];

x0 = log_beta0 - log(mu0);

ind_beta = 1;
ind_mu = 2;
ind_tau = 3;

funcs = set_funcs(...
    'sys_rhs', @(x, p) -p(ind_mu)*x(1,1)+exp(p(ind_beta)).*x(1,2).*exp(-x(1,2)), ...
    'sys_tau', @() ind_tau);

% choice of continuation parameter
ind_contpar = ind_beta;

%% Initialization of branch of non-trivial equilibria

maxstep=0.01; maxbound=10; 

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
    'contpar',[ind_mu,ind_beta],'dir',ind_mu,'step',0.01);
hbranch.method.continuation.plot=1;

figure(2);clf
hbranch=br_contn(funcs,hbranch,100);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,300);
xlabel('mu'); ylabel('log(beta)');

for jj = 1:size(hbranch.point,2)
    mu_hopf(jj) = hbranch.point(jj).parameter(ind_mu);
    log_beta_hopf(jj) = hbranch.point(jj).parameter(ind_beta);
end

figure(5); clf
plot(mu_hopf, exp(log_beta_hopf)./mu_hopf); hold on
axis([0 5 0 50])
xlabel('mu'); ylabel('beta/mu')
% savefig('db_hopf');



%% Branch off at  Hopf bifurcation
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');

ntst=80;

[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',1,'intervals',ntst,'degree',4,...
    'max_bound',[ind_contpar,maxbound],'max_step',[ind_contpar,maxstep],'radius',1,'step',1);
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
    'contpar',[ind_mu,ind_beta],'dir',ind_mu,'step',1e-2);
if ~suc
    error('fail',...
        ' initialization of period doubling failed');
end
figure(2);
pdbranch1=br_contn(pdfuncs,pdbranch1,100);
pdbranch=br_rvers(pdbranch1);
pdbranch=br_contn(pdfuncs,pdbranch,300);


for jj = 1:size(pdbranch.point,2)
    mu_pd(jj) = pdbranch.point(jj).parameter(ind_mu);
    log_beta_pd(jj) = pdbranch.point(jj).parameter(ind_beta);
end

figure(5); hold on
plot(mu_pd, exp(log_beta_pd)./mu_pd); hold on
axis([0 5 0 50])
xlabel('mu'); ylabel('beta/mu')
% savefig(['db_regions_ntst_',num2str(ntst)]);
% 
% save(['blowflies_DB_ntst_',num2str(ntst)])
