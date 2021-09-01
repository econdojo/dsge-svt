function par = ParTran(par,lb,ub,chi,con_unc)
% Function PARTRAN
%
% Purpose:    Perform parameter transformation
%
% Format:     par = ParTran(par,lb,ub,chi,con_unc)
%
% Input:      par       parameter value (N x dim)
%             lb        parameter lower bound (1 x dim)
%             ub        parameter upper bound (1 x dim)
%             chi       optimization tuning parameter
%             con_unc   constrained/unconstrained
%
% Output:     par       parameter value
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Parameter Transformation
%---------------------------------------------

% Initialization
dim = size(par,2);           % number of parameters

% Choose parameter transformation
if ~chi                      % no transformation (chi = 0)
    return
else
    switch con_unc
        case 'con'           % constrained
            for k = 1:dim
                if ~isinf(lb(k)) && isinf(ub(k))           % one-sided
                    par(:,k) = lb(k)+exp(chi*par(:,k));
                elseif ~isinf(lb(k)) && ~isinf(ub(k))      % two-sided
                    par(:,k) = lb(k)+(ub(k)-lb(k))./(1+exp(chi*par(:,k)));
                end
            end
        case 'unc'           % unconstrained
            for k = 1:dim
                if ~isinf(lb(k)) && isinf(ub(k))           % one-sided
                    par(:,k) = log(par(:,k)-lb(k))/chi;
                elseif ~isinf(lb(k)) && ~isinf(ub(k))      % two-sided
                    par(:,k) = log((ub(k)-lb(k))./(par(:,k)-lb(k))-1)/chi;
                end
            end
        otherwise
            error('Transformation does not exist!')
    end
end

%-------------------- END --------------------
