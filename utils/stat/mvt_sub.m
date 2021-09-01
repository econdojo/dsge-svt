function [mu1,Sigma1,v1] = mvt_sub(b1,b2,x2,mu,Sigma,v)
% Function MVT_SUB
%
% Purpose:    Compute conditional density of multivariate Student-t density
%             with arbitrarily ordered ordinates
%
% Format:     [mu1,Sigma1,v1] = mvt_sub(b1,b2,x2,mu,Sigma,v)
%
% Input:      b1        parameter index vector (1 x dim1)
%             b2        conditional parameter index vector (1 x dim2)
%             x2        conditional parameter row vector (1 x dim2)
%             mu        full mean vector (1 x dim)
%             Sigma     full scaling matrix (dim x dim); covariance = v/(v-2)*Sigma
%             v         full degrees of freedom
%
% Output:     mu1       conditional mean vector (1 x dim1)
%             Sigma1    conditional scaling matrix (dim1 x dim1)
%             v1        conditional degrees of freedom
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%        Conditional Student-t Density
%---------------------------------------------

% Initialization
b = union(b1,b2);
dim = length(mu);
if length(b)~=dim || max(b~=1:dim)
    error('Parameter index is incorrectly specified.')
else
    R = cholmod(Sigma); Sigma = R'*R;  % modified Cholesky factorization
end

% Check conditionality
if isempty(b2)
    % Unconditional density
    mu1 = mu(b1);
    Sigma1 = Sigma(b1,b1);
    v1 = v;
else
    % Conditional density
    Sigma11 = Sigma(b1,b1);
    Sigma12 = Sigma(b1,b2);
    Sigma22 = Sigma(b2,b2);
    mu1 = mu(b1)+(x2-mu(b2))/Sigma22*Sigma12';
    Sigma1 = Sigma11-Sigma12/Sigma22*Sigma12';

    % Choose density type
    if isinf(v)
        % Multivariate Normal
        v1 = Inf;
    else
        % Multivariate Student-t
        v1 = v+length(b2);
        Sigma1 = (v+(x2-mu(b2))/Sigma22*(x2-mu(b2))')/v1*Sigma1;
    end
end
%-------------------- END --------------------
