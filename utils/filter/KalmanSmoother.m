function [ss,Omega_ss] = KalmanSmoother(STE,fs,Omega_fs)
% Function KALMANSMOOTHER
%
% Purpose:    Compute smoothing densities; see
%             DeJong & Dave (2011), Structural Macroeconometrics
%
% Format:     [ss,Omega_ss] = KalmanSmoother(STE,fs,Omega_fs)
%
% Input:      STE       state transition equation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (diagonal)
%
%             fs        state filtering mean (from Kalman filter)
%             Omega_fs  state filtering covariance (from Kalman filter)
%
% Output:     ss(:,t)   E[s(t)|Y(:,1:T)] for t = 0,...,T
%             Omega_ss  Cov[s(t)|Y(:,1:T)] for t = 0,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%               Kalman Smoother
%---------------------------------------------

% Initialization
[dims,T] = size(fs);                             % numbers of state vars & periods
ss = [zeros(dims,T-1) fs(:,end)];                % state smoothing mean
Omega_ss = cat(3,zeros(dims,dims,T-1),Omega_fs(:,:,end));  % state smoothing covariance

% Transition covariance
if size(STE.Sigma_e,2)==1
    Sigma_e = repmat(STE.Sigma_e,1,T-1);         % constant
else
    Sigma_e = STE.Sigma_e;                       % time-varying
end

% Kalman smoother
for t = (T-1):-1:1
    % Period-t predictive density (degenerate)
    ps = STE.C+STE.G*fs(:,t);
    Omega_ps = STE.G*Omega_fs(:,:,t)*STE.G'+STE.M*diag(Sigma_e(:,t))*STE.M';
    
    % Period-t smoothing density
    %gain = (Omega_fs(:,:,t)*STE.G')*pinv(Omega_ps);   % Moore-Penrose pseudo-inverse
    gain = (Omega_fs(:,:,t)*STE.G')/Omega_ps;    % smoother gain
    ss(:,t) = fs(:,t)+gain*(ss(:,t+1)-ps);
    Omega_ss(:,:,t) = Omega_fs(:,:,t)+gain*(Omega_ss(:,:,t+1)-Omega_ps)*gain';
end

%-------------------- END --------------------
