function se = DistSmoother(Y,SSR,fs,Omega_fs) %#codegen
% Function DISTSMOOTHER
%
% Purpose:    Compute smoothed shocks; see
%             Koopman (1993), Disturbance smoother for state space models, Biometrika
%
% Format:     se = DistSmoother(Y,SSR,fs,Omega_fs)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (dime x T)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             fs        initial filtering mean
%             Omega_fs  initial filtering covariance
%
% Output:     se(:,t)   E[e(t)|Y(:,1:T)] for t = 1,...,T
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
% Updated: February 7, 2017

%% -------------------------------------------
%            Disturbance Smoother
%---------------------------------------------

% Initialization
T = size(Y,1);                              % number of periods
[dims,dime] = size(SSR.M);                  % numbers of states & shocks
dimy = length(SSR.D);                       % number of observables
pe = zeros(dimy,T);                         % prediction error
Omega_py = zeros(dimy,dimy,T);              % prediction covariance
gain = zeros(dims,dimy,T);                  % Kalman gain
se = zeros(dime,T);                         % smoothed shocks

% Kalman filter
for t = 1:T
    % Period-t state prediction
    ps = SSR.C+SSR.G*fs;
    Omega_ps = SSR.G*Omega_fs*SSR.G'+SSR.M*diag(SSR.Sigma_e(:,t))*SSR.M';
    
    % Period-t observable prediction
    py = SSR.D+SSR.Z*ps;
    pe(:,t) = Y(t,:)'-py;
    Omega_py(:,:,t) = SSR.Z*Omega_ps*SSR.Z'+diag(SSR.Sigma_u);
    
    % Period-t filtering density
    gain(:,:,t) = (Omega_ps*SSR.Z')/Omega_py(:,:,t);
    fs = ps+gain(:,:,t)*pe(:,t);
    Omega_fs = Omega_ps-gain(:,:,t)*SSR.Z*Omega_ps;
end

% Disturbance smoother
r = zeros(dims,1);
for t = T:-1:1
    r = SSR.Z'/Omega_py(:,:,t)*pe(:,t)+(SSR.G-SSR.G*gain(:,:,t)*SSR.Z)'*r;
    se(:,t) = diag(SSR.Sigma_e(:,t))*SSR.M'*r;
end

%-------------------- END --------------------
