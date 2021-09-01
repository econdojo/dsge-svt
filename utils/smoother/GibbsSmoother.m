function [fs,Omega_fs,loglik,S,R] = GibbsSmoother(Y,SSR,P,R,draw_S,draw_R)
% Function GIBBSSMOOTHER
%
% Purpose:    Sample states & regimes from smoothed distribtion; see
%             Carter & Kohn (1994), On Gibbs Sampling for State Space Models, Biometrika
%
% Format:     [fs,Omega_fs,loglik,S,R] = GibbsSmoother(Y,SSR,P,R,draw_S,draw_R)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = C(R_t) + G(R_t)*s(t-1) + e(t), V(e) = Sigma_e (R_t)
%                       y(t) = D(R_t) + Z(R_t)*s(t) + u(t), V(u) = Sigma_u (R_t)
%
%             P         transition matrix P(R_t=j|R_t-1=i) = P(i,j)
%             R(t)      current regime at time t
%             draw_S    update states based on regime draws
%             draw_R    update regimes based on state draws
%
% Output:     fs(:,t)   E[s(t)|Y(:,1:t)] for t = 0,...,T
%             Omega_fs  Cov[s(t)|Y(:,1:t)] for t = 0,...,T
%             loglik    log p(y(t)|Y(:,1:t-1)) for t = 1,...,T
%             S(:,t)    state draw at t = 1,...,T
%             R(t)      regime draw at t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: January 6, 2020

%% -------------------------------------------
%    Forward-Backward Simulation Smoother
%---------------------------------------------

% Initialization
T = size(Y,1);                         % number of periods
nr = size(P,1);                        % number of regimes
dims = size(SSR.G{1},1);               % number of states
fs = zeros(dims,T+1);                  % filtering mean
fs0 = zeros(dims,nr);
Omega_fs = zeros(dims,dims,T+1);       % filtering covariance
Omega_fs0 = zeros(dims,dims,nr);
fp = zeros(nr,T+1);                    % filtering probability
pp = zeros(nr,T);                      % predictive probability
S = zeros(dims,T+1);                   % state draws
loglik = zeros(T,1);                   % period log likelihood

% Initial distribution
for j = 1:nr
    fs0(:,j) = (eye(dims)-SSR.G{j})\SSR.C{j};
    Omega_fs0(:,:,j) = real(dlyap(SSR.G{j},SSR.Sigma_e{j}));
end

[V,D] = eig(P');
[~,id] = min(abs(diag(D)-1));
fp(:,1) = real(V(:,id));
fp(:,1) = fp(:,1)/sum(fp(:,1));

for j = 1:nr
    fs(:,1) = fs(:,1)+fp(j,1)*fs0(:,j);
end
    
for j = 1:nr
    Omega_fs(:,:,1) = Omega_fs(:,:,1)+fp(j,1)*(Omega_fs0(:,:,j)+(fs(:,1)-fs0(:,j))*(fs(:,1)-fs0(:,j))');
end

% Kalman filter
for t = 1:T
    % Period-(t-1) predictive density
    ps = SSR.C{R(t)}+SSR.G{R(t)}*fs(:,t);
    Omega_ps = SSR.G{R(t)}*Omega_fs(:,:,t)*SSR.G{R(t)}'+SSR.Sigma_e{R(t)};
    
    % Period-t log likelihood
    py = SSR.D{R(t)}+SSR.Z{R(t)}*ps;
    Omega_py = SSR.Z{R(t)}*Omega_ps*SSR.Z{R(t)}'+SSR.Sigma_u{R(t)};
    loglik(t) = mvt_pdf(Y(t,:),py',Omega_py,inf);
    
    % Period-t filtering density
    gain = (Omega_ps*SSR.Z{R(t)}')/Omega_py;     % Kalman gain
    fs(:,t+1) = ps+gain*(Y(t,:)'-py);
    Omega_fs(:,:,t+1) = Omega_ps-gain*SSR.Z{R(t)}*Omega_ps;
end

% Sample states
% Note: regard s(t+1) = C(R_t+1) + G(R_t+1)*s(t) + e(t+1) as measurement
% equation (error distribution must be non-degenerate)
if draw_S
    S(:,1) = fs(:,1);
    S(:,T+1) = mvt_rnd(fs(:,T+1)',Omega_fs(:,:,T+1),Inf,1)';
    for t = (T-1):-1:1
        % Period-t log likelihood
        py = SSR.C{R(t+1)}+SSR.G{R(t+1)}*fs(:,t+1);
        Omega_py = SSR.G{R(t+1)}*Omega_fs(:,:,t+1)*SSR.G{R(t+1)}'+SSR.Sigma_e{R(t+1)};

        % Period-t filtering density
        gain = (Omega_fs(:,:,t+1)*SSR.G{R(t+1)}')/Omega_py;    % Kalman gain
        fy = fs(:,t+1)+gain*(S(:,t+2)-py);
        Omega_fy = Omega_fs(:,:,t+1)-gain*SSR.G{R(t+1)}*Omega_fs(:,:,t+1);

        % Period-t state
        S(:,t+1) = mvt_rnd(fy',Omega_fy,Inf,1)';
    end
end

% Sample regimes
if draw_S && draw_R
    for t = 1:T
        % Period-(t-1) predictive probability
        pp(:,t) = P'*fp(:,t);

        % Period-t filtering probability
        for j = 1:nr
            fp(j,t+1) = exp(mvt_pdf(Y(t,:),(SSR.D{j}+SSR.Z{j}*S(:,t+1))',SSR.Sigma_u{j},Inf)...
                +mvt_pdf(S(:,t+1)',(SSR.C{j}+SSR.G{j}*S(:,t))',SSR.Sigma_e{j},Inf))*pp(j,t);
        end
        fp(:,t+1) = fp(:,t+1)/sum(fp(:,t+1));
    end

    R(T) = randsample(nr,1,true,fp(:,T+1));
    for t = (T-1):-1:1
        w = P(:,R(t+1)).*fp(:,t+1)/pp(R(t+1),t);
        R(t) = randsample(nr,1,true,w);
    end
elseif ~draw_S && draw_R
    error('Must sample states before regimes.')
end
S = S(:,2:end);

%-------------------- END --------------------
