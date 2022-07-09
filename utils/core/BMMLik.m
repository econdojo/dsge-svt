function loglik = BMMLik(SSR,nlag,mdof,para,data)
% Function BMMLIK
%
% Purpose:    Evaluate BMM log pseudo-likelihood
%
% Format:     loglik = BMMLik(SSR,nlag,mdof,para,data)
%
% Input:      SSR       state space representation (structure)
%             nlag      number of lags
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
[T,m] = size(data.M);
M = zeros(1,m);
n = length(E_y);
l = 1;
for k = 1:n
    if data.E_y(k)==0
        M(l) = E_y(k); l = l+1;
    elseif data.E_y(k)>0
        M(l) = E_y(k)+para(data.E_y(k)); l = l+1;
    end
end
for k = 1:nlag+1
    for i = 1:n
        for j = 1:n
            if data.E_yy(i,j,k)==0
                M(l) = E_yy(i,j,k); l = l+1;
            elseif data.E_yy(i,j,k)>0
                M(l) = E_yy(i,j,k)+para(data.E_yy(i,j,k)); l = l+1;
            end
        end
    end
end
V = NeweyWest_mex(data.M-repmat(M,T,1),-1);

% Evaluate Log Likelihood
if isinf(mdof)
    loglik = mvt_pdf_mex(mean(data.M),M,V,mdof);
else
    loglik = mvt_pdf_mex(mean(data.M),M,(mdof-2)/mdof*V,mdof);
end

%-------------------- END --------------------
