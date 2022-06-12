function loglik = BMMLik(SSR,nlag,mdof,para,data)
% Function BMMLIK
%
% Purpose:    Evaluate BMM log pseudo-likelihood
%
% Format:     loglik = VARLik(SSR,prior,nlag,data)
%
% Input:      SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e(t) (diagonal)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             nlag      number of VAR lags
%             mdof      kernel t degrees of freedom
%             para      vector of model parameters
%             data      structure with BMM data matrices
%
% Output:     loglik    BMM log pseudo-likelihood
%
% Written by Fei Tan, Saint Louis University
% Updated: June 10, 2022

%% -------------------------------------------
%           Evaluate Log Likelihood
%---------------------------------------------

% Assemble DSGE population moments
[E_y,E_yy] = Moment(SSR,nlag);
n = length(E_y);
M = [];
for k = 1:n
    if data.E_y(k)==0
        M = cat(2,M,E_y(k));
    elseif data.E_y(k)>0
        M = cat(2,M,E_y(k)+para(data.E_y(k)));
    end
end
for k = 1:nlag+1
    for i = 1:n
        for j = 1:n
            if data.E_yy(i,j,k)==0
                M = cat(2,M,E_yy(i,j,k));
            elseif data.E_yy(i,j,k)>0
                M = cat(2,M,E_yy(i,j,k)+para(data.E_yy(i,j,k)));
            end
        end
    end
end

% Evaluate Log Likelihood
if isinf(mdof)
    loglik = mvt_pdf_mex(data.M,M,data.V,mdof);
else
    loglik = mvt_pdf_mex(data.M,M,(mdof-2)/mdof*data.V,mdof);
end

%-------------------- END --------------------