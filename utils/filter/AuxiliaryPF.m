function [fs,loglik] = AuxiliaryPF(Y,SSR,s_old,w_par,sdof)
% Function AUXILIARYPF
%
% Purpose:    Compute filtered states & log likelihood; see
%             Herbst & Schorfheide (2015), Bayesian Estimation of DSGE Models
%
% Format:     [fs,loglik] = AuxiliaryPF(Y,SSR,s_old,w_par,sdof)
%
% Input:      Y(t,:)    observation at time t
%             SSR       state space representation (structure)
%
%                       s(t) = Transit(s(t-1),e(t)), E[ee'] = Sigma_e (diagonal)
%                       y(t) = Measure(s(t)) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             s_old     collection of initial particles (dims x N)
%             w_par     initial particle weights (N x 1)
%             sdof      shock t degrees of freedom
%
% Output:     fs(:,t)   E[s(t)|Y(:,1:t)]
%             loglik    log p(y(t)|Y(:,1:t-1)) for t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Auxiliary Particle Filter
%---------------------------------------------

% Initialization
resample = 0.50;                       % resample percentile
T = size(Y,1);                         % sample size
Y = [Y;Y(end,:)];                      % add dummy observation
[dims,N] = size(s_old);                % numbers of states, particles
fs = zeros(dims,T);                    % filtered states
loglik = zeros(T,1);                   % period log likelihood

% Sequential Monte Carlo
for t = 1:T
    % Period-t propagation
    [s_new,w_imp] = TransitPF(Y(t,:),SSR,s_old,sdof,'generic');
    
    % Period-t weights (two sets)
    w_adj = log(w_par)-Measure(Y(t,:),SSR,s_old,1);
    w_inc = w_adj+w_imp+Measure(Y(t,:),SSR,s_new,0);
    w_aux = w_inc+Measure(Y(t+1,:),SSR,s_new,1);
    
    % Period-t cumulative log likelihood
    num = log(sum(exp(w_inc)));
    den = log(sum(exp(w_adj)));
    loglik(t) = num-den;
    if isinf(loglik(t)) || isnan(loglik(t))
        loglik(t) = -Inf; fs = [];
        return
    end
    
    % Period-t filtering
    fs(:,t) = s_new*exp(w_inc)/exp(num);
    w_par = exp(w_aux)/sum(exp(w_aux));
    
    % Period-t resampling
    neff = 1/sum(w_par.^2);            % effective swarm size
    
    if neff < resample*N               % (adaptive) bootstrap particle filter
        index = Resample(w_par,'systematic');    % systematic
        s_old = s_new(:,index);        % update particles
        w_par = repmat(1/N,N,1);       % update weights
    else
        s_old = s_new;
    end
end

%-------------------- END --------------------
