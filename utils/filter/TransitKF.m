function [lamb,log_vol,w_imp] = TransitKF(y,SSR,fs_con,log_vol)
% Function TRANSITKF
%
% Purpose:    Generate new swarm from transition equation implied by linearized
%             DSGE model solution using approximately conditionally optimal
%             proposal distribution adapted to current observation via one-step
%             Kalman filter; see
%             Herbst & Schorfheide (2015), Bayesian Estimation of DSGE Models
%
% Format:     [lamb,log_vol,w_imp] = TransitKF(y,SSR,fs_con,log_vol)
%
% Input:      y         current observation (row vector)
%             SSR       state space representation (structure)
%             fs_con    conditional filtering mean (dims x N)
%             log_vol   log volatilities (dime x N)
%
% Output:     lamb      gamma precisions
%             log_vol   log volatilities
%             w_imp     log importance weights
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Particle Swarm Propagation
%---------------------------------------------

% Initialization
[dims,dime] = size(SSR.M);             % number of states & shocks
dimy = length(y);                      % number of observables
N = size(fs_con,2);                    % numbers of particles
w_imp = zeros(N,1);                    % log importance weights

% Generate new swarm of volatilities (not adapted; naive forward propagation)
if SSR.sv
    shock = mvt_rnd(zeros(1,dime),diag(SSR.Sigma_v),Inf,N)';
    log_vol = repmat(SSR.A,1,N)+SSR.B*log_vol+shock;
end

% Generate new swarm of gamma precisions (adapted; generic importance sampler)
if ~isinf(SSR.sdof)
    % Construct expanded linear SSR
    C = [SSR.C;zeros(dime,1)];
    G = [SSR.G zeros(dims,dime);zeros(dime,dims+dime)];
    M = [SSR.M;eye(dime)];
    D = SSR.D;
    Z = [SSR.Z zeros(dimy,dime)];

    % Period-(t-1) predictive density
    ps = repmat(C,1,N)+G*[fs_con;zeros(dime,N)];
    Omega_ps = M*diag(SSR.Sigma_e(:,1))*M';

    % Period-t log likelihood
    py = repmat(D,1,N)+Z*ps;
    Omega_py = Z*Omega_ps*Z'+diag(SSR.Sigma_u);

    % Period-t filtering density
    fs = ps+(Omega_ps*Z')/Omega_py*(repmat(y',1,N)-py);
    %Omega_fs = Omega_ps-(Omega_ps*Z')/Omega_py*Z*Omega_ps;

    % Generate new swarm of gamma precisions
    fs = fs((dims+1):end,:);       % extract filtered/smoothed shocks
    lamb = gamrnd((SSR.sdof+dime)/2,2./(sum((fs'/diag(SSR.Sigma_e(:,1))).*fs',2)+SSR.sdof))';

    % Compute importance weights
    w_imp = log(gampdf(lamb,SSR.sdof/2,2/SSR.sdof))...
        -log(gampdf(lamb,(SSR.sdof+dime)/2,2./(sum((fs'/diag(SSR.Sigma_e(:,1))).*fs',2)+SSR.sdof)'));
else
    lamb = ones(1,N);
end

%-------------------- END --------------------
