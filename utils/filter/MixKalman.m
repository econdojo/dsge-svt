function [fs,Omega_fs,loglik] = MixKalman(Y,SSR,fs_init,Omega_fs_init,N) %#codegen
% Function MIXKALMAN
%
% Purpose:    Compute filtered states & log likelihood; see
%             Chen & Liu (2000), Mixture Kalman filters
%
% Format:     [fs,Omega_fs,loglik] = MixKalman(Y,SSR,fs_init,Omega_fs_init,N)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (dime x T)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%                       h(t) = A + B*h(t-1) + v(t), E[vv'] = Sigma_v (diagonal)
%                       sdof = shock t degrees of freedom
%                       sv   = false if using given volatility Sigma_e(t):
%                            - time-varying if Gaussian shock
%                            - time-invariant if t shock
%
%             fs_init        initial filtering mean
%             Omega_fs_init  initial filtering covariance
%             N         number of particles
%
% Output:     fs(:,t)   E[s(t)|Y(:,1:t)] for t = 0,...,T
%             Omega_fs  Cov[s(t)|Y(:,1:t)] for t = 0,...,T
%             loglik    log p(y(t)|Y(:,1:t-1)) for t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%            Mixture Kalman Filter
%---------------------------------------------

% Choose filter
if isinf(SSR.sdof) && ~SSR.sv
    % Conventional Kalman filter
    [fs,Omega_fs,loglik] = KalmanFilter(Y,SSR,fs_init,Omega_fs_init);
    return
else
    % Initialization
    resample = 0.50;                        % resample percentile
    T = size(Y,1);                          % sample size
    [dims,dime] = size(SSR.M);              % number of states & shocks
    fs = [fs_init zeros(dims,T)];           % filtered states
    fs_con = repmat(fs_init,1,N);           % conditional filtering mean
    Omega_fs = cat(3,Omega_fs_init,zeros(dims,dims,T));    % filtered covariance
    Omega_fs_con = repmat(Omega_fs_init,1,1,N);  % conditional filtering covariance
    log_vol = repmat((eye(dime)-SSR.B)\SSR.A,1,N);    % initial log volatilities
    w_par = repmat(1/N,N,1);                % initial particle weights
    w_inc = zeros(N,1);                     % initial incremental weights
    loglik = zeros(T,1);                    % period log likelihood
    if SSR.sv
        SSR.Sigma_e = repmat(exp(log_vol(:,1)),1,T);
    end
end

% Mixture Kalman filter
for t = 1:T
    % Period-t propagation
    [lamb,log_vol,w_imp] = TransitKF(Y(t,:),SSR,fs_con,log_vol);
    if SSR.sv
        Sigma_e_con = exp(log_vol)./repmat(lamb,dime,1);
    else
        Sigma_e_con = repmat(SSR.Sigma_e(:,1),1,N)./repmat(lamb,dime,1);
    end
    
    % Period-t conditional Kalman
    for k = 1:N
        % Period-(t-1) predictive density
        ps = SSR.C+SSR.G*fs_con(:,k);
        Omega_ps = SSR.G*Omega_fs_con(:,:,k)*SSR.G'+SSR.M*diag(Sigma_e_con(:,k))*SSR.M';

        % Period-t log likelihood
        py = SSR.D+SSR.Z*ps;
        Omega_py = SSR.Z*Omega_ps*SSR.Z'+diag(SSR.Sigma_u);

        % Cumulative log likelihood
        w_inc(k) = log(w_par(k))+w_imp(k)+mvt_pdf(Y(t,:),py',Omega_py,inf);
        
        % Period-t filtering density
        gain = (Omega_ps*SSR.Z')/Omega_py;  % Kalman gain
        fs_con(:,k) = ps+gain*(Y(t,:)'-py);
        Omega_fs_con(:,:,k) = Omega_ps-gain*SSR.Z*Omega_ps;
    end
    
    % Period-t cumulative log likelihood
    loglik(t) = log(sum(exp(w_inc)));
    if isinf(loglik(t)) || isnan(loglik(t))
        loglik(t) = -Inf; fs = []; Omega_fs = [];
        return
    end
    
    % Period-t filtering
    w_par = exp(w_inc-loglik(t));
    
    % Period-t resampling
    neff = 1/sum(w_par.^2);            % effective swarm size
    if neff < resample*N               % (adaptive) bootstrap particle filter
        index = Resample(w_par,'systematic');    % systematic
        fs_con = fs_con(:,index);      % update conditional densities
        Omega_fs_con = Omega_fs_con(:,:,index);
        log_vol = log_vol(:,index);    % update log volatilities
        w_par = repmat(1/N,N,1);       % update weights
    end
    
    % Filtered states & covariances
    fs(:,t+1) = fs_con*w_par;
    %Omega_fs(:,:,t+1) = sum(bsxfun(@times,Omega_fs_con,permute(w_old,[3 2 1])),3); % costly
end

%-------------------- END --------------------
