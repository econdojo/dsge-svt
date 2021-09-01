function loglik = LogLik(SSR,prior,nlag,data)
% Function LOGLIK
%
% Purpose:    Evaluate DSGE log marginal likelihood; see
%             Del Negro & Schorfheide (2004), Priors from general equilibrium
%             models for VARs, International Economic Review
%
% Format:     loglik = LogLik(SSR,prior,nlag,data)
%
% Input:      SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (diagonal)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             prior     DSGE prior weight
%             nlag      number of VAR lags
%             data      structure with VAR data matrices
%
% Output:     loglik    DSGE log marginal likelihood
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%           Evaluate Log Likelihood
%---------------------------------------------

% Initialization
[T,n] = size(data.Y);             % time span & number of data variables
k = 1+n*nlag;                     % number of regressors

% Assemble DSGE population moments
[E_y,E_yy_k] = Moment(SSR,nlag);
E_yx = [E_y reshape(E_yy_k(:,:,2:end),n,n*nlag)];
E_xx = zeros(n*nlag);
for i = 1:nlag
    for j = i:nlag
        E_xx(((i-1)*n+1):(i*n),((j-1)*n+1):(j*n)) = E_yy_k(:,:,j-i+1);
    end
end
E_xx = triu(E_xx)+triu(E_xx,1)';
E_xx = [1 repmat(E_y',1,nlag);repmat(E_y,nlag,1) E_xx];

% Check positive definiteness
[~,indef,~] = cholmod_mex(E_xx);
if indef
    loglik = -Inf;
    return
else
    % VAR prior MLE
    Phi_prior = E_xx\E_yx';
    E_uu_prior = E_yy_k(:,:,1)-E_yx*Phi_prior;
    
    % DSGE log marginal likelihood
    if isinf(prior)               % VAR approx of DSGE
        U = data.Y-data.X*Phi_prior;
        loglik = -n*T/2*log(2*pi)-T/2*logdet(E_uu_prior,'chol')-trace(E_uu_prior\U'*U)/2;
    else
        % VAR posterior MLE
        Phi_post = (prior*T*E_xx+data.X'*data.X)\(prior*T*E_yx'+data.X'*data.Y);
        E_uu_post = (prior*T*E_yy_k(:,:,1)+data.Y'*data.Y-(prior*T*E_yx+data.Y'*data.X)*Phi_post)/((prior+1)*T);
        
        ld1 = logdet(prior*T*E_xx+data.X'*data.X,'chol');
        ld2 = logdet((prior+1)*T*E_uu_post,'chol');
        ld3 = logdet(prior*T*E_xx,'chol');
        ld4 = logdet(prior*T*E_uu_prior,'chol');
        loglik = -n/2*ld1-((prior+1)*T-k)/2*ld2+n/2*ld3+(prior*T-k)/2*ld4-n*T/2*log(2*pi)+n*T/2*log(2)...
            +sum(gammaln(((prior+1)*T-k+1-(1:n))/2))-sum(gammaln((prior*T-k+1-(1:n))/2));
    end
end

%-------------------- END --------------------
