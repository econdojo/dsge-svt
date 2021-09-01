function [fs,loglik] = ParticleFilter(Y,SSR,swarm,w_par,sdof)
% Function PARTICLEFILTER
%
% Purpose:    Compute filtered states & log likelihood; see
%             DeJong & Dave (2011), Structural Macroeconometrics
%
% Format:     [fs,loglik] = ParticleFilter(Y,SSR,swarm,w_par,sdof)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = Transit(s(t-1),e(t)), E[ee'] = Sigma_e
%                       y(t) = Measure(s(t)) + u(t), E[uu'] = Sigma_u
%
%             swarm     collection of initial particles (dims x N)
%             w_par     initial particle weights (N x 1)
%             sdof      shock t degrees of freedom
%
% Output:     fs(:,t)   E[s(t)|Y(:,1:t)]
%             loglik    log p(y(t)|Y(:,1:t-1)) for t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%           Generic Particle Filter
%---------------------------------------------

% Initialization
resample = 0.50;                       % resample percentile
T = size(Y,1);                         % sample size
[dims,N] = size(swarm);                % numbers of states, particles
fs = zeros(dims,T);                    % filtered states
loglik = zeros(T,1);                   % period log likelihood

% Particle filter
for t = 1:T
    % Period-t propagation
    [swarm,w_imp] = TransitPF(Y(t,:),SSR,swarm,sdof,'generic');
    
    % Period-t cumulative log likelihood
    w_inc = log(w_par)+w_imp+Measure(Y(t,:),SSR,swarm,0);
    loglik(t) = log(sum(exp(w_inc)));
    if isinf(loglik(t)) || isnan(loglik(t))
        loglik(t) = -Inf; fs = [];
        return
    end
    
    % Period-t filtering
    w_par = exp(w_inc-loglik(t));
    
    % Period-t resampling
    neff = 1/sum(w_par.^2);            % effective swarm size
    
    if neff < resample*N               % (adaptive) bootstrap particle filter
        index = Resample(w_par,'systematic');    % systematic
        swarm = swarm(:,index);        % update particles
        w_par = repmat(1/N,N,1);       % update weights
    end
    
    fs(:,t) = swarm*w_par;             % filtered states
end

%-------------------- END --------------------
