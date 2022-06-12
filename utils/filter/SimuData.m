function [Y,shock] = SimuData(SSR,T,state,log_vol) %#codegen
% Function SIMUDATA
%
% Purpose:    Simulate data from state space representation
%
% Format:     [Y,shock] = SimuData(SSR,T,state,log_vol)
%
% Input:      SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (dime x T)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%                       h(t) = A + B*h(t-1) + v(t), E[vv'] = Sigma_v (diagonal)
%                       sdof = shock t degrees of freedom
%                       sv   = false if using given volatility Sigma_e(t):
%                            - time-varying if Gaussian shock
%                            - time-invariant if t shock
%
%             T         number of periods
%             state     initial state
%             log_vol   initial log volatility
%
% Output:     Y         data row vectors (T x dimy)
%             shock     shock row vectors (T x dime)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%               Data Simulation
%---------------------------------------------

% Initialization
[dims,dime] = size(SSR.M);             % numbers of states & shocks
dimy = length(SSR.D);                  % number of observables
Y = zeros(T,dimy);                     % observable data
shock = zeros(T,dime);                 % structural shock
if nargin<3
    state = (eye(dims)-SSR.G)\SSR.C;   % initial state
end
if nargin<4
    log_vol = (eye(dime)-SSR.B)\SSR.A; % initial log volatility
end

% Transition covariance
if SSR.sv
    SSR.Sigma_e = zeros(dime,T);
end

% Measurement covariance
if any(SSR.Sigma_u)
    err_u = mvt_rnd(zeros(1,dimy),diag(SSR.Sigma_u),Inf,T);
else
    err_u = zeros(T,dimy);             % no measurement error
end

% Volatility covariance
if any(SSR.Sigma_v)
    err_v = mvt_rnd(zeros(1,dime),diag(SSR.Sigma_v),Inf,T);
else
    err_v = zeros(T,dime);             % constant volatility
end

% Simulate data
for t = 1:T
    if SSR.sv
        log_vol = SSR.A+SSR.B*log_vol+err_v(t,:)';
        SSR.Sigma_e(:,t) = exp(log_vol);
    end
    shock(t,:) = mvt_rnd(zeros(1,dime),diag(SSR.Sigma_e(:,t)),SSR.sdof,1);
    state = SSR.C+SSR.G*state+SSR.M*shock(t,:)';
    Y(t,:) = SSR.D'+state'*SSR.Z'+err_u(t,:);
end

%-------------------- END --------------------
