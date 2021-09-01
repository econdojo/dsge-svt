function [fs,Omega_fs,loglik] = KalmanFilter(Y,SSR,fs_init,Omega_fs_init)
% Function KALMANFILTER
%
% Purpose:    Compute filtering densities & log likelihood; see
%             DeJong & Dave (2011), Structural Macroeconometrics
%
% Format:     [fs,Omega_fs,loglik] = KalmanFilter(Y,SSR,fs_init,Omega_fs_init)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (dime x T)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             fs_init        initial filtering mean
%             Omega_fs_init  initial filtering covariance
%
% Output:     fs(:,t)   E[s(t)|Y(:,1:t)] for t = 0,...,T
%             Omega_fs  Cov[s(t)|Y(:,1:t)] for t = 0,...,T
%             loglik    log p(y(t)|Y(:,1:t-1)) for t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%                Kalman Filter
%---------------------------------------------

% Initialization
T = size(Y,1);                              % number of periods
dims = length(fs_init);                     % number of state variables
fs = [fs_init zeros(dims,T)];               % filtering mean
Omega_fs = cat(3,Omega_fs_init,zeros(dims,dims,T));   % filtering covariance
loglik = zeros(T,1);                        % period log likelihood

% Kalman filter
for t = 1:T
    % Period-(t-1) predictive density
    ps = SSR.C+SSR.G*fs(:,t);
    Omega_ps = SSR.G*Omega_fs(:,:,t)*SSR.G'+SSR.M*diag(SSR.Sigma_e(:,t))*SSR.M';
    
    % Period-t log likelihood
    py = SSR.D+SSR.Z*ps;
    Omega_py = SSR.Z*Omega_ps*SSR.Z'+diag(SSR.Sigma_u);
    
    % Cumulative log likelihood
    loglik(t) = mvt_pdf(Y(t,:),py',Omega_py,inf);
    
    % Period-t filtering density
    gain = (Omega_ps*SSR.Z')/Omega_py;      % Kalman gain
    fs(:,t+1) = ps+gain*(Y(t,:)'-py);
    Omega_fs(:,:,t+1) = Omega_ps-gain*SSR.Z*Omega_ps;
end

%-------------------- END --------------------
