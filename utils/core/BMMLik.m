function loglik = BMMLik(E_y,E_yy,mdof,para,data)
% Function BMMLIK
%
% Purpose:    Evaluate BMM log pseudo-likelihood
%
% Format:     loglik = VARLik(SSR,prior,nlag,data)
%
% Input:      E_y       1st-order moments: E[y(t)] = E_y
%             E_yy      2nd-order raw moments: E[y(t)y(t-k)'] = E_yy(:,:,k+1)
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
M = zeros(1,length(data.M));
n = length(E_y);
nlag = size(E_yy,3)-1;
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

% Evaluate Log Likelihood
if isinf(mdof)
    loglik = mvt_pdf(data.M,M,data.V,mdof);
else
    loglik = mvt_pdf(data.M,M,(mdof-2)/mdof*data.V,mdof);
end

%-------------------- END --------------------