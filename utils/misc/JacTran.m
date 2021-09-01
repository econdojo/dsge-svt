function Jac = JacTran(par,lb,ub,chi)
% Function JACTRAN
%
% Purpose:    Compute Jacobian of parameter transformation
%
% Format:     Jac = JacTran(par,lb,ub,chi)
%
% Input:      par       unconstrained parameter value (N x dim)
%             lb        parameter lower bound (1 x dim)
%             ub        parameter upper bound (1 x dim)
%             chi       optimization tuning parameter
%
% Output:     Jac       Jacobian diagonal (dim x N)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Jacobian of Transformation
%---------------------------------------------

% Initialization
[N,dim] = size(par);         % numbers of copies & parameters
Jac = ones(dim,N);           % Jacobian diagonals (column)

% Choose parameter transformation
if ~chi                      % no transformation (chi = 0)
    return
else
    for k = 1:dim
        if ~isinf(lb(k)) && isinf(ub(k))         % one-sided
            Jac(k,:) = chi*exp(chi*par(:,k)');
        elseif ~isinf(lb(k)) && ~isinf(ub(k))    % two-sided
            Jac(k,:) = -(ub(k)-lb(k))*chi*exp(chi*par(:,k)')./(1+exp(chi*par(:,k)')).^2;
        end
    end
end

%-------------------- END --------------------
