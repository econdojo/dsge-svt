function test_bench
% Function TEST_BENCH
%
% Purpose:    DSGE model analysis & generate MEX functions
%
% Format:     test_bench
%
% Input:      none
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: January 4, 2018

clear
close all
clc

%% -------------------------------------------
%             DSGE Model Analysis
%---------------------------------------------

% Set path
modpath = ['user' filesep 'as07'];
datpath = ['user' filesep 'as07' filesep 'data.txt'];
savepath = ['user' filesep 'as07'];
addpath(genpath('utils'),genpath('mex'),modpath);

% State space representation (SSR)
%load([savepath filesep 'chain_init.mat'],'xopt','pid')% load post mode
%load([savepath filesep 'tarb_full.mat'],'stat')  % load post mean
[P,V] = ParVar;                   % load parameters & variables
Y = importdata(datpath);          % import data
%P.mod.para(pid) = xopt';
%P.mod.para(P.mod.svp) = stat(P.mod.svp,1);
para = P.mod.para;
G0 = zeros(V.mod.nvar);
G1 = zeros(V.mod.nvar);
Psi = zeros(V.mod.nvar,V.mod.nshock);
Pi = zeros(V.mod.nvar,V.mod.nfore);
CC = zeros(V.mod.nvar,1);
Sigma_u = zeros(V.mod.ndata,1);
Z = zeros(V.mod.ndata,V.mod.nvar);
D = zeros(V.mod.ndata,1);
ssp = zeros(P.mod.nssp,1); %#ok<NASGU>
j = 0; %#ok<NASGU>
user_ssp
user_mod
SSR.Sigma_u = Sigma_u;
SSR.Z = Z;
SSR.D = D;
[SSR.G,SSR.C,SSR.M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);
fprintf('Existence  =  %d,  Uniqueness  =  %d\n',eu(1),eu(2));

% Stochastic volatility
T = 200;                    % number of periods
SSR.sdof = 5;                   % Student-t degrees of freedom
SSR.sv = ~isempty(P.mod.svp);     % stochastic volatility
if isempty(P.mod.svp)
    Sigma_e = para(P.mod.svp_ss).^2;
    SSR.Sigma_v = zeros(V.mod.nshock,1);
    SSR.B = zeros(V.mod.nshock);
    SSR.A = zeros(V.mod.nshock,1);
else
    Sigma_e = exp(para(P.mod.svp_ss));
    SSR.Sigma_v = para(P.mod.svp_vv);
    SSR.B = diag(para(P.mod.svp_ar));
    SSR.A = (eye(V.mod.nshock)-SSR.B)*para(P.mod.svp_ss);
end
SSR.Sigma_e = repmat(Sigma_e,1,T);

% Plot impulse response
% plotvar = {'y', 'pi', 'R'};       % variables
% irfplot(SSR,plotvar,V,80);        % execute plot

% Simulate data
rng(2311);                        % repeatable random number generation
[Y,~] = SimuData(SSR,T);
dlmwrite([savepath filesep 'simu_data.txt'],Y,'delimiter','\t','precision',10);

%% -------------------------------------------
%           Likelihood Evaluation
%---------------------------------------------

% Auxiliary particle filter (only t shock)
N = 5000*((V.mod.nvar>20)+1);
swarm = zeros(V.mod.nvar,N);
SSR.Sigma_e = Sigma_e;
shock = mvt_rnd_mex(zeros(1,V.mod.nshock),diag(Sigma_e),SSR.sdof,N)';
swarm = repmat(SSR.C,1,N)+SSR.G*swarm+SSR.M*shock;
y_prob = Measure(Y(1,:),SSR,swarm,1);
w_par = mean(y_prob);
w_par = exp(y_prob-w_par)/sum(exp(y_prob-w_par));
tic;
[~,loglik_apf] = AuxiliaryPF(Y,SSR,swarm,w_par,SSR.sdof);
fprintf('Log likelihood by auxiliary particle filter  =  %.3f\n',sum(loglik_apf));
toc;

% Mixture Kalman filter (t shock & stochastic vol)
SSR.Sigma_e = repmat(Sigma_e,1,T);
fs = (eye(V.mod.nvar)-SSR.G)\SSR.C;
Omega_fs = real(dlyap(SSR.G,SSR.M*diag(Sigma_e)*SSR.M'));
N = 1000*V.mod.nshock;            % number of particles
tic;
[~,~,loglik_mkf] = MixKalman(Y,SSR,fs,Omega_fs,N);
fprintf('Log likelihood by mixture Kalman filter  =  %.3f\n',sum(loglik_mkf));
toc;

% Kalman filter & smoother (Gaussian shock & constant vol)
SSR.C = [SSR.C;zeros(V.mod.nshock,1)]; % expanded linear SSR
SSR.G = [SSR.G zeros(V.mod.nvar,V.mod.nshock);zeros(V.mod.nshock,V.mod.nvar+V.mod.nshock)];
SSR.M = [SSR.M;eye(V.mod.nshock)];
SSR.Z = [SSR.Z zeros(V.mod.ndata,V.mod.nshock)];
fs = (eye(V.mod.nvar+V.mod.nshock)-SSR.G)\SSR.C;
Omega_fs = real(dlyap(SSR.G,SSR.M*diag(Sigma_e)*SSR.M'));
tic;
[~,~,loglik_kf] = KalmanFilter(Y,SSR,fs,Omega_fs);
fprintf('Log likelihood by standard Kalman filter  =  %.3f\n',sum(loglik_kf));
toc;
DistSmoother(Y,SSR,fs,Omega_fs);

%-------------------- END --------------------
